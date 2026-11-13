"""Numbat (Gao et al. 2023): haplotype-aware CNVs from scRNA-seq, or scMultiome with binned
RNA + ATAC counts (Numbat-multiome: the "genes" are genomic bins named ``chr1:0-1042457``).

Reads a Numbat ``out_dir`` (or a folder holding ``numbat_out/``), at its final iteration ``i``
(the highest with ``clone_post_i.tsv``, ``segs_consensus_i.tsv`` and ``geno_i.tsv``):

- ``clone_post_i.tsv``: per cell ``clone_opt``, ``compartment_opt`` (tumor / normal), ``p_cnv``
- ``segs_consensus_i.tsv``: the CNV segments (``seg_cons``, ``seg_start``, ``seg_end``,
  ``cnv_state``: neu, del, amp, loh, bamp, bdel)
- ``geno_i.tsv``: cells x segments, the posterior probability that the cell carries the CNV
- ``gexp_roll_wide.tsv.gz``: cells x genes (or bins), the smoothed expression log ratio

The bins are the columns of ``gexp_roll_wide`` when they are genomic bins; for a gene-level run,
1 Mb windows. Layers:

- ``p_cnv``: on each bin, the probability of the most likely CNV segment over it (0 outside)
- ``total_cn`` / ``major_cn`` / ``minor_cn``: where ``p_cnv`` >= ``min_p`` (0.9), that
  segment's state as allele copies, otherwise neutral 1:1 (del 1:0, amp 2:1, loh 2:0, bamp 2:2,
  bdel 0:0): relative states, the way Numbat's heatmap shows them, not absolute copy numbers
- ``logr``: the smoothed expression ratio (bin-level runs only: a gene-level run's outputs do
  not hold the genes' positions)

Cells: ``clone``, ``compartment``, ``p_cnv_cell`` (Numbat's p_cnv), ``p_opt``, ``cnv_burden``.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from ..model import CNVDataset, Modality, Semantics, chr_sort_key
from . import Reader, register

STATES = {
    "neu": (1, 1),
    "del": (1, 0),
    "amp": (2, 1),
    "loh": (2, 0),
    "bamp": (2, 2),
    "bdel": (0, 0),
}
BIN_NAME = re.compile(r"^(chr)?([0-9XYM]+):(\d+)-(\d+)$")
WINDOW = 1_000_000
MAX_BINS = 8000
#: a segment's state is painted on a cell when Numbat is this sure the cell carries it. Per-cell
#: segment posteriors are often near 0.5 (little evidence per cell, e.g. on short segments):
#: Numbat's own heatmap fades them out with p_cnv; 0.9 keeps the confident ones.
MIN_P = 0.9


def _final(d: Path) -> int | None:
    its = []
    for p in d.glob("clone_post_*.tsv"):
        m = re.match(r"clone_post_(\d+)\.tsv$", p.name)
        if (
            m
            and (d / f"segs_consensus_{m[1]}.tsv").is_file()
            and (d / f"geno_{m[1]}.tsv").is_file()
        ):
            its.append(int(m[1]))
    return max(its) if its else None


def runs_under(path: Path) -> list[Path]:
    """Numbat output folders at `path` or up to two levels below it."""
    found = []
    for d in [path, *path.glob("*"), *path.glob("*/*")]:
        if d.is_dir() and not d.name.startswith(".") and _final(d) is not None:
            found.append(d)
    return found


def _out_dir(path: Path) -> Path:
    runs = runs_under(path)
    if not runs:
        raise ValueError(
            f"no Numbat output (clone_post_N.tsv, segs_consensus_N.tsv, geno_N.tsv) in {path}"
        )
    if len(runs) > 1:
        names = ", ".join(str(r.relative_to(path)) for r in runs)
        raise ValueError(
            f"{path} holds several Numbat runs ({names}); point --input at one of them"
        )
    return runs[0]


def _detect(path: Path) -> bool:
    return (
        path.is_dir()
        and _final(path) is not None
        or (path / "numbat_out").is_dir()
        and (_final(path / "numbat_out") is not None)
    )


def _bin_grid(names: list[str]) -> pd.DataFrame | None:
    """Bins from gexp_roll_wide's columns when they are genomic bins (Numbat-multiome)."""
    rows = []
    for n in names:
        m = BIN_NAME.match(n)
        if not m:
            return None
        rows.append(("chr" + m[2], int(m[3]), int(m[4]), n))
    return pd.DataFrame(rows, columns=["chr", "start", "end", "name"])


def _coarsen(bins: pd.DataFrame, X: np.ndarray | None) -> tuple[pd.DataFrame, np.ndarray | None]:
    """Merge neighbouring bins of a chromosome until there are at most MAX_BINS (the heatmap's
    width); X (cells x bins) is averaged over the merged bins."""
    k = int(np.ceil(len(bins) / MAX_BINS))
    if k <= 1:
        return bins, X
    grp = bins.groupby("chr", sort=False).cumcount() // k
    key = bins["chr"] + ":" + grp.astype(str)
    codes, uniq = pd.factorize(key)
    nb = bins.groupby(codes).agg(chr=("chr", "first"), start=("start", "min"), end=("end", "max"))
    if X is not None:
        S = np.zeros((X.shape[0], len(uniq)))
        C = np.zeros((X.shape[0], len(uniq)))
        fin = np.isfinite(X)
        np.add.at(S.T, codes, np.where(fin, X, 0).T)
        np.add.at(C.T, codes, fin.T)
        with np.errstate(invalid="ignore", divide="ignore"):
            X = np.where(C > 0, S / C, np.nan)
    return nb.reset_index(drop=True), X


def read_numbat(
    path: Path, genome: str = "hg38", modality: str | None = None, min_p: float = MIN_P
) -> CNVDataset:
    d = _out_dir(path)
    i = _final(d)
    cp = pd.read_csv(d / f"clone_post_{i}.tsv", sep="\t", dtype={"cell": str})
    cp = cp.drop_duplicates("cell").set_index("cell")
    cells = cp.index.tolist()
    segs = pd.read_csv(d / f"segs_consensus_{i}.tsv", sep="\t", dtype={"CHROM": str})
    if "seg_cons" not in segs:
        segs["seg_cons"] = segs["seg"]
    segs = segs.drop_duplicates("seg_cons").set_index("seg_cons")
    geno = pd.read_csv(d / f"geno_{i}.tsv", sep="\t", dtype={"cell": str}).set_index("cell")
    geno = geno.reindex(index=cells)
    geno = geno.loc[:, geno.columns.isin(segs.index)]

    # the bins, and the smoothed expression on them
    gexp_file = d / "gexp_roll_wide.tsv.gz"
    bins, logr, binned = None, None, False
    if gexp_file.is_file():
        head = pd.read_csv(gexp_file, sep="\t", nrows=0).columns.tolist()
        cell_col = head[0]
        bins = _bin_grid(head[1:])
        if bins is not None:  # gene-level runs: gene positions are not in Numbat's outputs
            binned = True
            types = {c: np.float32 for c in head[1:]} | {cell_col: str}
            g = pd.read_csv(gexp_file, sep="\t", dtype=types, engine="c")
            g = g.set_index(cell_col).reindex(index=cells)
            order = np.lexsort((bins["start"].to_numpy(), bins["chr"].map(chr_sort_key).to_numpy()))
            bins = bins.iloc[order].reset_index(drop=True)
            logr = g[bins["name"]].to_numpy(dtype=np.float32)
            del g
            bins = bins.drop(columns="name")
    if bins is None:
        chroms = segs.assign(c="chr" + segs["CHROM"].astype(str)).groupby("c")["seg_end"].max()
        rows = [
            (c, k * WINDOW, (k + 1) * WINDOW)
            for c in sorted(chroms.index, key=chr_sort_key)
            for k in range(int(np.ceil(chroms[c] / WINDOW)))
        ]
        bins = pd.DataFrame(rows, columns=["chr", "start", "end"])
    bins, logr = _coarsen(bins, logr)

    # each segment's probability on the bins it overlaps; a bin takes its most likely segment
    n, nb = len(cells), len(bins)
    P = np.zeros((n, nb), dtype=np.float32)
    S = np.full((n, nb), -1, dtype=np.int8)  # index into state_names
    state_names = list(STATES)
    chr_b = bins["chr"].to_numpy()
    s_b, e_b = bins["start"].to_numpy(), bins["end"].to_numpy()
    G = geno.to_numpy(dtype=np.float32)
    for j, seg in enumerate(geno.columns):
        r = segs.loc[seg]
        st = str(r.get("cnv_state_post", r["cnv_state"]))
        if st not in STATES or st == "neu":
            continue
        m = (chr_b == "chr" + str(r["CHROM"])) & (s_b < r["seg_end"]) & (e_b > r["seg_start"])
        if not m.any():
            continue
        p = np.nan_to_num(G[:, j], nan=0.0)[:, None]
        better = p > P[:, m]
        P[:, m] = np.where(better, p, P[:, m])
        S[:, m] = np.where(better, state_names.index(st), S[:, m])
    missing = geno.isna().all(axis=1).to_numpy()  # cells Numbat gave no genotype
    P[missing] = np.nan

    A = np.ones((n, nb), dtype=np.float32)
    B = np.ones((n, nb), dtype=np.float32)
    call = (P >= min_p) & (S >= 0)
    for k, st in enumerate(state_names):
        m = call & (S == k)
        A[m], B[m] = STATES[st]
    A[missing] = B[missing] = np.nan

    obs = pd.DataFrame(index=cells)
    if "clone_opt" in cp:
        obs["clone"] = pd.Categorical(cp["clone_opt"].astype("string"))
    if "compartment_opt" in cp:
        obs["compartment"] = pd.Categorical(cp["compartment_opt"].astype("string"))
    for src, dst in (("p_cnv", "p_cnv_cell"), ("p_opt", "p_opt")):
        if src in cp:
            obs[dst] = pd.to_numeric(cp[src], errors="coerce").to_numpy()
    with np.errstate(invalid="ignore"):
        obs["cnv_burden"] = np.where(missing, np.nan, call.mean(axis=1))

    layers = {
        "total_cn": A + B,
        "major_cn": np.fmax(A, B),
        "minor_cn": np.fmin(A, B),
        "p_cnv": P,
    }
    sem = {
        "total_cn": Semantics.INTEGER_CN,
        "major_cn": Semantics.ALLELE_CN,
        "minor_cn": Semantics.ALLELE_CN,
        "p_cnv": Semantics.PROBABILITY,
    }
    if logr is not None:
        layers["logr"], sem["logr"] = logr, Semantics.CONTINUOUS_RATIO
    mod = Modality(modality) if modality else (Modality.MULTIOME if binned else Modality.SCRNA)
    # cells ordered by their CNV posteriors, as Numbat builds its phylogeny from geno (the
    # painted states are too sparse to order the cells by); both views use it
    dendros = {}
    if 3 <= n <= 20000 and not missing.all():
        from ..analysis import ward_manhattan

        d_ = ward_manhattan(np.where(missing[:, None], 0.0, P))
        method = "ward.D2, Manhattan on the CNV posteriors (p_cnv)"
        dendros = {
            "total": d_ | {"layers": ["p_cnv"], "method": method},
            "allele": d_ | {"layers": ["p_cnv"], "method": method},
        }
    root = d.parent if d.name == "numbat_out" else d
    return CNVDataset.build(
        cells=obs,
        bins=bins,
        layers={k: v.astype(np.float32) for k, v in layers.items()},
        semantics=sem,
        modality=mod,
        tool="Numbat",
        genome=genome,
        default_layer="total_cn",
        extra_uns={
            "source": root.name,
            "numbat_iteration": i,
            "numbat_binned": binned,
            "numbat_min_p": min_p,
            "dendrograms": dendros,
            # expression-based calls: the DNA QC metrics (MAPD, residuals, ...) do not apply
            "qc": {"version": 1, "computed": ["cnv_burden"]},
        },
    )


def _read(path: Path, **opts) -> CNVDataset:
    return read_numbat(path, **opts)


register(
    Reader("numbat", "scRNA", _detect, _read, "Numbat out_dir (clone_post_N.tsv, geno_N.tsv, ...)")
)
