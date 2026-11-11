"""Chromosome sizes, bin grids, and projecting segment calls onto bins."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .model import chr_sort_key

CHROM_SIZES: dict[str, dict[str, int]] = {
    "hg38": {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    },
    "hg19": {
        "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
        "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
        "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
        "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
        "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
        "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566,
    },
}  # fmt: skip


def normalise_chr(name: object, style: str = "chr") -> str:
    """'1', 'chr1', 'Chr1' → 'chr1' (style='chr') or '1' (style='bare')."""
    bare = str(name).removeprefix("chr").removeprefix("Chr").removeprefix("CHR")
    if bare in ("23",):
        bare = "X"
    if bare in ("24",):
        bare = "Y"
    return f"chr{bare}" if style == "chr" else bare


def uniform_bins(
    genome: str = "hg38", size: int = 1_000_000, chromosomes: list[str] | None = None
) -> pd.DataFrame:
    """A fixed-width grid over the genome: columns chr, start (0-based), end (exclusive)."""
    sizes = CHROM_SIZES[genome]
    chroms = chromosomes or [c for c in sizes if c != "chrY"]
    rows = []
    for c in sorted(chroms, key=chr_sort_key):
        n = sizes[c]
        starts = np.arange(0, n, size, dtype=np.int64)
        ends = np.minimum(starts + size, n)
        rows.append(pd.DataFrame({"chr": c, "start": starts, "end": ends}))
    return pd.concat(rows, ignore_index=True)


def segments_to_bins(
    segments: pd.DataFrame,
    bins: pd.DataFrame,
    value_cols: list[str],
    cell_col: str = "cell",
    cells: list[str] | None = None,
) -> dict[str, np.ndarray]:
    """Project per-cell segments onto a bin grid.

    ``segments`` has columns chr, start, end, ``cell_col`` and the ``value_cols``. Each bin takes
    the value of the segment covering its midpoint (NaN where none does). Returns one
    cells × bins float32 matrix per value column, rows in the order of ``cells``.
    """
    cells = cells if cells is not None else list(pd.unique(segments[cell_col].astype(str)))
    row_of = {c: i for i, c in enumerate(cells)}
    out = {v: np.full((len(cells), len(bins)), np.nan, dtype=np.float32) for v in value_cols}
    mids = ((bins["start"].to_numpy() + bins["end"].to_numpy()) // 2).astype(np.int64)
    bin_chr = bins["chr"].astype(str).to_numpy()
    chr_bins = {c: np.flatnonzero(bin_chr == c) for c in pd.unique(bin_chr)}

    seg = segments.copy()
    seg["_chr"] = seg["chr"].map(
        lambda x: normalise_chr(x, "chr" if bin_chr[0].startswith("chr") else "bare")
    )
    seg["_row"] = seg[cell_col].astype(str).map(row_of)
    seg = seg.dropna(subset=["_row"])
    for (c, row), g in seg.groupby(["_chr", "_row"], sort=False):
        idx = chr_bins.get(c)
        if idx is None or not len(idx):
            continue
        m = mids[idx]
        g = g.sort_values("start")
        s = g["start"].to_numpy()
        e = g["end"].to_numpy()
        # segment whose start is the last one <= midpoint, then check it reaches the midpoint
        k = np.searchsorted(s, m, side="right") - 1
        ok = (k >= 0) & (e[np.clip(k, 0, None)] >= m)
        for v in value_cols:
            vals = g[v].to_numpy(dtype=np.float32)
            out[v][int(row), idx[ok]] = vals[k[ok]]
    return out
