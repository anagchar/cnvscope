"""scphasedcn (SPLONGGET): allele-specific copy number from single-cell long reads.

Reads a caller folder written by ``write_scphasedcn()`` / ``call_scn.R`` (or a sample folder
holding ``caller/``):

- ``cells.tsv``: the called cells (cell type, depth stratum, null cell, multiplet; tumour mode:
  clone, posterior, doublet_like, ...)
- normal mode, ``calls.tsv``: called cell x region rows (arm, chromosome, 5-50 Mb tiles,
  events) with a state (loss 1:0, gain 2:1, cnloh 2:0, gain_as 3:0, amp, del 0:0)
- tumour mode, ``clone_segments.tsv``: allele-specific segments (``state`` = "nA:nB") per clone

Layers on 1 Mb windows: ``total_cn``, ``major_cn`` / ``minor_cn``. Normal mode paints each
cell's reportable calls over the neutral state (broad regions first, nested ones over them);
tumour mode paints each assigned cell with its clone's segments (unassigned cells stay empty).

When the extraction folder is found (``extract/`` beside the caller folder, or the ``dir`` in
the run notes), its 100 kb bins add ``read_counts`` and ``hap_a_frac`` (haplotype-A fraction
of the phased molecules, a BAF) summed to the same windows.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..model import CNVDataset, Modality, Semantics, chr_sort_key
from . import Reader, register

WINDOW = 1_000_000
CALLER_FILES = ("cells.tsv", "calls.tsv", "thresholds.tsv")
#: normal-mode state -> (nA, nB) on a diploid chromosome (caller_states(), cn_grid())
STATES = {
    "neutral": (1, 1),
    "loss": (1, 0),
    "gain": (2, 1),
    "cnloh": (2, 0),
    "gain_as": (3, 0),
    "amp": (3, 1),
    "del": (0, 0),
}
#: on a haploid chromosome (a male's X / Y): state_label()
HAPLOID = {"neutral": (1, 0), "loss": (0, 0), "del": (0, 0), "gain": (2, 0), "amp": (3, 0)}
#: broad regions first, so nested calls paint over them
CLASS_ORDER = {"chrom": 0, "arm": 1, "50Mb": 2, "20Mb": 3, "10Mb": 4, "5Mb": 5, "event": 6}


def _caller_dir(path: Path) -> Path | None:
    for d in (path, path / "caller"):
        if d.is_dir() and all((d / f).is_file() for f in CALLER_FILES):
            if (d / "params.yaml").is_file() or (d / "result.rds").is_file():
                return d
    return None


def _detect(path: Path) -> bool:
    d = _caller_dir(path)
    if d is None:
        return False
    with open(d / "calls.tsv") as f:
        head = f.readline().rstrip("\n").split("\t")
    return {"cell", "region", "class", "state", "level"} <= set(head)


def _notes(d: Path) -> dict[str, str]:
    """`- key: value` lines of NOTES.md (inputs of the run)."""
    p = d / "NOTES.md"
    if not p.is_file():
        return {}
    out = {}
    for line in p.read_text(errors="replace").splitlines():
        m = re.match(r"^\s*-\s*([A-Za-z_]+):\s*(.+?)\s*$", line)
        if m and m.group(1) not in out:
            out[m.group(1)] = m.group(2)
    return out


def _extract_dir(caller: Path, notes: dict[str, str]) -> Path | None:
    cands = []
    if "dir" in notes:
        p = Path(notes["dir"])
        cands += [p] if p.is_absolute() else [caller.parent / p, caller.parent.parent / p]
    cands += [caller.parent / "extract", caller / "extract"]
    for c in cands:
        if c.is_dir() and any(c.glob("bins_*.h5ad")):
            return c
    return None


def _windows(chroms: dict[str, int]) -> pd.DataFrame:
    rows = []
    for c in sorted(chroms, key=chr_sort_key):
        n = int(np.ceil(chroms[c] / WINDOW))
        for k in range(n):
            rows.append((c, k * WINDOW, min((k + 1) * WINDOW, chroms[c])))
    return pd.DataFrame(rows, columns=["chr", "start", "end"])


def _haploid(chrom: str, sex: str | None) -> bool:
    return sex == "male" and chrom.removeprefix("chr") in ("X", "Y")


def _paint(seg: pd.DataFrame, bins: pd.DataFrame, A: np.ndarray, B: np.ndarray, row: int) -> None:
    """Write segments (chrom, start, end, nA, nB) onto row `row`: windows whose midpoint falls
    inside a segment take its state."""
    mid = (bins["start"].to_numpy() + bins["end"].to_numpy()) / 2
    chr_ = bins["chr"].to_numpy()
    for c, s, e, a, b in seg[["chrom", "start", "end", "nA", "nB"]].itertuples(index=False):
        m = (chr_ == c) & (mid >= s) & (mid < e)
        A[row, m], B[row, m] = a, b


def read_scphasedcn(path: Path, genome: str = "hg38") -> CNVDataset:
    d = _caller_dir(path)
    if d is None:
        raise ValueError(f"no scphasedcn caller folder (cells.tsv, calls.tsv, ...) in {path}")
    notes = _notes(d)
    cells_t = pd.read_csv(d / "cells.tsv", sep="\t", dtype={"cell": str}).set_index("cell")
    cells = cells_t.index.astype(str).tolist()
    seg_file = d / "clone_segments.tsv"
    tumour = seg_file.is_file() and seg_file.stat().st_size > 0 and "clone" in cells_t
    sex = notes.get("sex")

    # the genome: chromosome extents from the scored regions (or the extraction's bins)
    ext = _extract_dir(d, notes)
    reg = pd.read_csv(d / "regions.tsv", sep="\t") if (d / "regions.tsv").is_file() else None
    chroms: dict[str, int] = {}
    if reg is not None and len(reg):
        chroms = reg.groupby("chrom")["end"].max().astype(np.int64).to_dict()
    if tumour:
        seg = pd.read_csv(seg_file, sep="\t", dtype={"clone": str})
        for c, e in seg.groupby("chrom")["end"].max().items():
            chroms[c] = max(chroms.get(c, 0), int(e))
    if not chroms:
        raise ValueError(f"{d}: no regions.tsv to take the chromosomes from")
    bins = _windows(chroms)
    nb, n = len(bins), len(cells)
    hap = np.array([_haploid(c, sex) for c in bins["chr"]])
    A = np.full((n, nb), np.nan, dtype=np.float32)
    B = np.full((n, nb), np.nan, dtype=np.float32)

    if tumour:
        seg[["nA", "nB"]] = seg["state"].str.split(":", expand=True).astype(int)
        by_clone = {k: g for k, g in seg.groupby("clone")}
        assigned = cells_t["assigned"] if "assigned" in cells_t else pd.Series(True, index=cells)
        for i, c in enumerate(cells):
            cl = cells_t.at[c, "clone"]
            if pd.isna(cl) or not bool(assigned.get(c, True)) or str(cl) not in by_clone:
                continue
            _paint(by_clone[str(cl)], bins, A, B, i)
        mode = "tumour"
    else:
        A[:] = np.where(hap, 1, 1)
        B[:] = np.where(hap, 0, 1)
        calls = pd.read_csv(d / "calls.tsv", sep="\t", dtype={"cell": str})
        keep = calls["level"].eq("cell") & calls["called"].astype(str).eq("True")
        if "reportable" in calls:
            keep &= calls["reportable"].astype(str).eq("True")
        calls = calls[keep & calls["cell"].isin(cells_t.index)].copy()
        calls["_o"] = calls["class"].map(CLASS_ORDER).fillna(9)
        calls = calls.sort_values(["_o", "end"], ascending=[True, True])
        row = {c: i for i, c in enumerate(cells)}
        for c, g in calls.groupby("cell", sort=False):
            g = g.sort_values("_o", kind="stable")
            ab = [
                (HAPLOID if _haploid(ch, sex) else STATES).get(st, (np.nan, np.nan))
                for ch, st in zip(g["chrom"], g["state"], strict=True)
            ]
            g = g.assign(nA=[x[0] for x in ab], nB=[x[1] for x in ab]).dropna(subset=["nA"])
            _paint(g, bins, A, B, row[c])
        mode = "normal"

    layers = {
        "total_cn": A + B,
        "major_cn": np.fmax(A, B),
        "minor_cn": np.fmin(A, B),
    }
    sem = {
        "total_cn": Semantics.INTEGER_CN,
        "major_cn": Semantics.ALLELE_CN,
        "minor_cn": Semantics.ALLELE_CN,
    }
    if ext is not None:
        extra = _depth(ext, cells, bins)
        layers |= extra
        sem |= {k: (Semantics.COUNT if k == "read_counts" else Semantics.BAF) for k in extra}

    obs = cells_t.copy()
    obs = obs.loc[:, obs.notna().any()]  # columns the run left empty
    for c in obs.columns:
        if obs[c].dtype == object:
            obs[c] = pd.Categorical(obs[c].astype("string"))
    return CNVDataset.build(
        cells=obs,
        bins=bins,
        layers={k: v.astype(np.float32) for k, v in layers.items()},
        semantics=sem,
        modality=Modality.SCDNA,
        tool="scphasedcn",
        genome=genome,
        default_layer="total_cn",
        extra_uns={
            "source": (d.parent.name if d.name == "caller" else d.name),
            "scphasedcn_mode": mode,
            "sex": sex or "",
            "window_bp": WINDOW,
            "extraction": str(ext) if ext else "",
        },
    )


def _depth(ext: Path, cells: list[str], bins: pd.DataFrame) -> dict[str, np.ndarray]:
    """Read counts and the haplotype-A fraction of `ext`/bins_*.h5ad, summed to `bins`."""
    import anndata as ad
    import scipy.sparse as sp

    f = sorted(ext.glob("bins_*.h5ad"))[0]
    a = ad.read_h5ad(f)
    idx = pd.Index(a.obs_names.astype(str))
    rows = idx.get_indexer(cells)
    have = rows >= 0
    v = a.var
    key = pd.MultiIndex.from_arrays([bins["chr"].to_numpy(), bins["start"].to_numpy() // WINDOW])
    src = pd.MultiIndex.from_arrays(
        [v["chrom"].astype(str).to_numpy(), (v["start"].to_numpy() // WINDOW).astype(np.int64)]
    )
    col = key.get_indexer(src)  # window of each 100 kb bin (-1: outside the grid)
    ok = col >= 0
    agg = sp.csr_matrix(
        (np.ones(ok.sum()), (np.flatnonzero(ok), col[ok])), shape=(a.n_vars, len(bins))
    )

    def summed(m) -> np.ndarray:
        out = np.full((len(cells), len(bins)), np.nan, dtype=np.float64)
        m = sp.csr_matrix(m) if not sp.issparse(m) else m
        out[have] = np.asarray((m[rows[have]] @ agg).todense())
        return out

    out = {"read_counts": summed(a.X)}
    if "hapA" in a.layers and "hapB" in a.layers:
        ha, hb = summed(a.layers["hapA"]), summed(a.layers["hapB"])
        with np.errstate(invalid="ignore", divide="ignore"):
            out["hap_a_frac"] = np.where(ha + hb > 0, ha / (ha + hb), np.nan)
    return out


def _read(path: Path, **opts) -> CNVDataset:
    return read_scphasedcn(path, **opts)


register(
    Reader(
        "scphasedcn",
        "scDNA",
        _detect,
        _read,
        "scphasedcn caller folder (cells.tsv, calls.tsv / clone_segments.tsv)",
    )
)
