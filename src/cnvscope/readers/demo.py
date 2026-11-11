"""A synthetic scDNA dataset with clones, allele-specific copy number and noisy logR.

Used for development and tests until real tool fixtures are in place.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..genome import uniform_bins
from ..model import CNVDataset, Modality, Semantics
from . import Reader, register

# (chr, start fraction, end fraction, major, minor) events per clone, on top of 1|1
CLONES = {
    "A": [("chr3", 0.5, 1.0, 2, 1), ("chr8", 0.0, 1.0, 2, 1), ("chr17", 0.0, 0.5, 1, 0)],
    "B": [("chr3", 0.5, 1.0, 2, 1), ("chr10", 0.0, 1.0, 1, 0), ("chr20", 0.0, 1.0, 3, 1)],
    "C": [("chr7", 0.0, 1.0, 2, 1), ("chr9", 0.2, 0.6, 0, 0), ("chr13", 0.0, 1.0, 2, 0)],
}


def make_demo(n_cells: int = 300, bin_size: int = 1_000_000, seed: int = 0) -> CNVDataset:
    rng = np.random.default_rng(seed)
    bins = uniform_bins("hg38", bin_size)
    nb = len(bins)
    clone_names = ["normal", *CLONES]
    weights = np.array([0.15, 0.35, 0.3, 0.2])
    clones = rng.choice(clone_names, size=n_cells, p=weights)

    major = np.ones((n_cells, nb), dtype=np.float32)
    minor = np.ones((n_cells, nb), dtype=np.float32)
    chr_arr = bins["chr"].to_numpy()
    for clone, events in CLONES.items():
        rows = np.flatnonzero(clones == clone)
        for c, f0, f1, ma, mi in events:
            idx = np.flatnonzero(chr_arr == c)
            a, b = idx[0] + int(f0 * len(idx)), idx[0] + int(f1 * len(idx))
            major[np.ix_(rows, np.arange(a, b))] = ma
            minor[np.ix_(rows, np.arange(a, b))] = mi
    # a few private events per cell
    for i in range(n_cells):
        for _ in range(rng.poisson(0.6)):
            a = rng.integers(0, nb - 20)
            ln = rng.integers(3, 20)
            major[i, a : a + ln] = rng.choice([1, 2, 3])
            minor[i, a : a + ln] = rng.choice([0, 1])
    total = major + minor
    noise = rng.normal(0, rng.uniform(0.08, 0.25, size=(n_cells, 1)), size=(n_cells, nb))
    ploidy = total.mean(axis=1, keepdims=True)
    logr = (np.log2(np.maximum(total, 0.05) / ploidy) + noise).astype(np.float32)
    # some missing bins
    logr[rng.random(logr.shape) < 0.002] = np.nan

    depth = rng.uniform(15, 60, size=(n_cells, 1))
    reads = rng.poisson(depth * np.maximum(total, 0.1) / ploidy).astype(np.float32)
    ci_major = rng.uniform(0.1, 0.45, size=(n_cells, 1)).astype(np.float32)
    ci_minor = rng.uniform(0.1, 0.35, size=(n_cells, 1)).astype(np.float32)
    cell_type = np.where(
        clones == "normal",
        rng.choice(["Immune", "Stromal", "Endothelial"], size=n_cells),
        "Epithelial",
    )
    cells = pd.DataFrame(
        {
            "clone": pd.Categorical(clones),
            "cell_type": pd.Categorical(cell_type),
            "coverage": (depth[:, 0] / 8).round(1),
        },
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )
    ds = CNVDataset.build(
        cells=cells,
        bins=bins,
        layers={
            "total_cn": total,
            "major_cn": major,
            "minor_cn": minor,
            "logr": logr,
            "read_counts": reads,
            "major_cn_lower": np.maximum(0, major - ci_major),
            "major_cn_upper": major + ci_major,
            "minor_cn_lower": np.maximum(0, minor - ci_minor),
            "minor_cn_upper": minor + ci_minor,
        },
        semantics={
            "total_cn": Semantics.INTEGER_CN,
            "major_cn": Semantics.ALLELE_CN,
            "minor_cn": Semantics.ALLELE_CN,
            "logr": Semantics.CONTINUOUS_RATIO,
            "read_counts": Semantics.COUNT,
            "major_cn_lower": Semantics.ALLELE_CN,
            "major_cn_upper": Semantics.ALLELE_CN,
            "minor_cn_lower": Semantics.ALLELE_CN,
            "minor_cn_upper": Semantics.ALLELE_CN,
        },
        modality=Modality.SCDNA,
        tool="demo",
        genome="hg38",
        default_layer="total_cn",
    )
    from ..prepare import prepare

    return prepare(ds)


def _detect(path: Path) -> bool:
    return str(path) == "demo" and not path.exists()


def _read(path: Path, **opts) -> CNVDataset:
    return make_demo(**opts)


register(Reader("demo", "scDNA", _detect, _read, "Synthetic scDNA dataset"))
