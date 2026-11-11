"""Everything a dataset gets after it is read and before it is shown."""

from __future__ import annotations

from .analysis import cluster_cells
from .model import CNVDataset
from .qc import compute_qc


def prepare(ds: CNVDataset, cluster: bool = True) -> CNVDataset:
    """Add missing dendrograms (total and allele-specific) and QC metrics."""
    if cluster:
        cluster_cells(ds)
    compute_qc(ds)
    return ds
