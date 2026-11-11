"""CHISEL (Zaccaria & Raphael 2021): allele- and haplotype-specific copy number per cell.

Reads the output directory of ``chisel`` (or ``calls/calls.tsv`` itself):

- ``calls/calls.tsv``: one row per bin and cell, columns ``#CHR START END CELL NORM_COUNT
  COUNT RDR A_COUNT B_COUNT BAF CLUSTER HAP_CN CORRECTED_HAP_CN``. ``HAP_CN`` is ``a|b``,
  the copy numbers of the two haplotypes; ``CORRECTED_HAP_CN`` is used when present.
- ``clones/mapping.tsv`` (optional): ``#CELL CLUSTER CLONE``; ``None`` = no clone.

Layers: ``total_cn`` (a + b), ``major_cn``/``minor_cn`` (max/min of the two haplotypes),
``hap_a_cn``/``hap_b_cn`` (the haplotypes as phased), ``rdr``, ``baf``, ``read_counts``
(``COUNT``), and ``cn_scale``: RDR on the copy-number scale, RDR × the cell's mean copy
number / its mean RDR, so bins can be drawn around the calls.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..model import CNVDataset, Modality, Semantics, chr_sort_key
from . import Reader, register


def _calls_path(path: Path) -> Path | None:
    if path.is_file() and path.name.endswith(".tsv"):
        return path
    for p in (path / "calls" / "calls.tsv", path / "calls.tsv"):
        if p.is_file():
            return p
    return None


def _detect(path: Path) -> bool:
    p = _calls_path(path)
    if p is None:
        return False
    with open(p) as f:
        head = f.readline()
    return head.startswith("#CHR") and "HAP_CN" in head and "CELL" in head


def _matrix(rows: np.ndarray, cols: np.ndarray, vals: np.ndarray, shape) -> np.ndarray:
    m = np.full(shape, np.nan, dtype=np.float32)
    m[rows, cols] = vals
    return m


def read_chisel(path: Path, genome: str = "hg38") -> CNVDataset:
    calls = _calls_path(path)
    if calls is None:
        raise ValueError(f"no CHISEL calls.tsv under {path}")
    df = pd.read_csv(calls, sep="\t", dtype={"#CHR": str, "CELL": str})
    df = df.rename(columns={"#CHR": "chr"})
    hap_col = "CORRECTED_HAP_CN" if "CORRECTED_HAP_CN" in df else "HAP_CN"
    ab = df[hap_col].astype(str).str.split("|", n=1, expand=True)
    a = pd.to_numeric(ab[0], errors="coerce").to_numpy(np.float32)
    b = pd.to_numeric(ab[1], errors="coerce").to_numpy(np.float32)

    bins = df[["chr", "START", "END"]].drop_duplicates()
    bins = bins.assign(_k=bins["chr"].map(chr_sort_key)).sort_values(["_k", "START"])
    bins = bins.drop(columns="_k").reset_index(drop=True)
    bin_index = pd.MultiIndex.from_frame(bins)
    col = bin_index.get_indexer(pd.MultiIndex.from_frame(df[["chr", "START", "END"]]))
    cells = list(pd.unique(df["CELL"]))
    row = pd.Index(cells).get_indexer(df["CELL"])
    shape = (len(cells), len(bins))

    hap_a, hap_b = _matrix(row, col, a, shape), _matrix(row, col, b, shape)
    total = hap_a + hap_b
    rdr = _matrix(row, col, df["RDR"].to_numpy(np.float32), shape)
    with np.errstate(invalid="ignore", divide="ignore"):
        scale = np.nanmean(total, axis=1, keepdims=True) / np.nanmean(rdr, axis=1, keepdims=True)
    layers = {
        "total_cn": total,
        "major_cn": np.fmax(hap_a, hap_b),
        "minor_cn": np.fmin(hap_a, hap_b),
        "hap_a_cn": hap_a,
        "hap_b_cn": hap_b,
        "rdr": rdr,
        "cn_scale": (rdr * scale).astype(np.float32),
        "baf": _matrix(row, col, df["BAF"].to_numpy(np.float32), shape),
        "read_counts": _matrix(row, col, df["COUNT"].to_numpy(np.float32), shape),
    }
    sem = {
        "total_cn": Semantics.INTEGER_CN,
        "major_cn": Semantics.ALLELE_CN,
        "minor_cn": Semantics.ALLELE_CN,
        "hap_a_cn": Semantics.ALLELE_CN,
        "hap_b_cn": Semantics.ALLELE_CN,
        "rdr": Semantics.RATIO,
        "cn_scale": Semantics.CONTINUOUS_CN,
        "baf": Semantics.BAF,
        "read_counts": Semantics.COUNT,
    }

    first = df.drop_duplicates("CELL").set_index("CELL")
    obs = pd.DataFrame(index=cells)
    obs["cluster"] = pd.Categorical(first.loc[cells, "CLUSTER"].astype(str))
    root = calls.parent.parent if calls.parent.name == "calls" else calls.parent
    mapping = root / "clones" / "mapping.tsv"
    if mapping.is_file():
        mp = pd.read_csv(mapping, sep="\t", dtype=str).rename(columns={"#CELL": "CELL"})
        clone = mp.set_index("CELL")["CLONE"].reindex(cells)
        clone = clone.where(~clone.isin(["None", "NA", ""]))
        obs["clone"] = pd.Categorical(clone)

    return CNVDataset.build(
        cells=obs,
        bins=bins.rename(columns={"START": "start", "END": "end"}),
        layers=layers,
        semantics=sem,
        modality=Modality.SCDNA,
        tool="CHISEL",
        genome=genome,
        default_layer="total_cn",
        extra_uns={"source": root.name, "hap_cn_column": hap_col},
    )


def _read(path: Path, **opts) -> CNVDataset:
    return read_chisel(path, **opts)


register(Reader("chisel", "scDNA", _detect, _read, "CHISEL output directory (calls/calls.tsv)"))
