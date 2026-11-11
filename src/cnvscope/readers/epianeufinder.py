"""epiAneufinder (Ramakrishnan et al. 2023): gains and losses from scATAC-seq read depth.

Reads the run's output folder (``epiAneufinder_results/``, or the folder above it):

- ``results_table.tsv``: R ``write.table`` output (space-separated, with row names), columns
  ``seq start end cell-<barcode> ...``, one row per genomic window, values 0 = loss, 1 =
  normal, 2 = gain.

Layer ``state`` (-1 loss, 0 normal, +1 gain). Cells (barcodes without the ``cell-`` prefix):
``cnv_burden`` (fraction of windows not normal, epiAneufinder's ``cnv_burden_per_cell``),
``gain_fraction`` and ``loss_fraction``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..model import CNVDataset, Modality, Semantics, chr_sort_key
from . import Reader, register

TABLE = "results_table.tsv"
PREFIX = "cell-"


def _table(path: Path) -> Path | None:
    if path.is_file():
        return path if path.name.endswith(".tsv") else None
    for p in (path / TABLE, path / "epiAneufinder_results" / TABLE):
        if p.is_file():
            return p
    return None


def _detect(path: Path) -> bool:
    t = _table(path)
    if t is None:
        return False
    with open(t) as f:
        head = f.readline().split()
    return head[:3] == ["seq", "start", "end"] and any(h.startswith(PREFIX) for h in head[3:9])


def read_epianeufinder(path: Path, genome: str = "hg38") -> CNVDataset:
    t = _table(path)
    if t is None:
        raise ValueError(f"no epiAneufinder {TABLE} in {path}")
    with open(t) as f:
        line = f.readline()
        head, first = line.split(), f.readline().split()
    sep = "\t" if "\t" in line else " "
    cols = [h for h in head[3:] if h.startswith(PREFIX)]
    # R's write.table adds row names: then data rows have one more field than the header
    names = ["_row", *head] if len(first) == len(head) + 1 else head
    raw = pd.read_csv(
        t,
        sep=sep,
        skiprows=1,
        header=None,
        names=names,
        dtype={"seq": str},
        na_values=["NA"],
        engine="c",
    )
    bins = pd.DataFrame(
        {
            "chr": raw["seq"].astype(str),
            "start": raw["start"].astype(np.int64),
            "end": raw["end"].astype(np.int64),
        }
    )
    order = np.lexsort((bins["start"].to_numpy(), bins["chr"].map(chr_sort_key).to_numpy()))
    bins = bins.iloc[order].reset_index(drop=True)
    calls = raw[cols].to_numpy(dtype=np.float32)[order].T  # cells x windows
    state = calls - 1  # 0/1/2 -> -1/0/+1
    fin = np.isfinite(state)
    n = np.maximum(fin.sum(axis=1), 1)
    cells = [c[len(PREFIX) :] for c in cols]
    obs = pd.DataFrame(
        {
            "cnv_burden": (np.where(fin, state != 0, False).sum(axis=1) / n),
            "gain_fraction": (np.where(fin, state > 0, False).sum(axis=1) / n),
            "loss_fraction": (np.where(fin, state < 0, False).sum(axis=1) / n),
        },
        index=cells,
    )
    root = t.parent.parent if t.parent.name == "epiAneufinder_results" else t.parent
    return CNVDataset.build(
        cells=obs,
        bins=bins,
        layers={"state": state.astype(np.float32)},
        semantics={"state": Semantics.TERNARY_STATE},
        modality=Modality.SCATAC,
        tool="epiAneufinder",
        genome=genome,
        default_layer="state",
        extra_uns={"source": root.name},
    )


def _read(path: Path, **opts) -> CNVDataset:
    return read_epianeufinder(path, **opts)


register(
    Reader(
        "epianeufinder",
        "scATAC",
        _detect,
        _read,
        "epiAneufinder output folder (epiAneufinder_results/results_table.tsv)",
    )
)
