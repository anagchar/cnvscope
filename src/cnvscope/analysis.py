"""Clustering of cells, computed from the canonical dataset.

Matches the ASCAT.sc Explorer and ASCAT.sc's R plots: Ward clustering (R's ``ward.D2``) of
the Manhattan distances between profiles, missing values set to -1. scipy's ``ward`` on a
precomputed distance vector applies the same Lance-Williams update as ``ward.D2``.
"""

from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from .model import CNVDataset, Semantics


def scipy_to_hclust(Z: np.ndarray) -> dict:
    """Convert a scipy linkage matrix to R's hclust ``merge``/``height``.

    hclust: row i merges two items; negative = leaf (-1-based leaf index), positive = the
    1-based row of an earlier merge. scipy: indices < n are leaves, n + i is merge row i.
    """
    n = Z.shape[0] + 1
    merge = np.empty((n - 1, 2), dtype=np.int64)
    for i in range(n - 1):
        for j in range(2):
            a = int(Z[i, j])
            merge[i, j] = -(a + 1) if a < n else (a - n + 1)
    return {"merge": merge, "height": Z[:, 2].astype(np.float64)}


def hclust_order(merge: np.ndarray) -> np.ndarray:
    """Leaf order of an hclust merge matrix (left child first), as row indices."""
    merge = np.asarray(merge, dtype=np.int64)
    if not len(merge):
        return np.zeros(1, dtype=np.int64)
    out: list[int] = []
    stack = [len(merge)]
    while stack:
        code = stack.pop()
        if code < 0:
            out.append(-code - 1)
            continue
        left, right = merge[code - 1]
        stack.append(int(right))
        stack.append(int(left))
    return np.asarray(out, dtype=np.int64)


def ward_manhattan(x: np.ndarray) -> dict:
    """Ward.D2 on Manhattan distances of the rows of ``x`` (NaN → -1)."""
    x = np.where(np.isnan(x), -1.0, x.astype(np.float64))
    Z = linkage(pdist(x, metric="cityblock"), method="ward")
    h = scipy_to_hclust(Z)
    h["order"] = leaves_list(Z).astype(np.int64)
    return h


def cluster_cells(ds: CNVDataset, max_cells: int = 20000, force: bool = False) -> CNVDataset:
    """Add the ``total`` dendrogram (default layer) and, when the dataset has major and minor
    copy number, the ``allele`` dendrogram (major and minor stacked). Existing ones are kept
    unless ``force``. Skipped for fewer than 3 or more than ``max_cells`` cells.
    """
    n = ds.adata.n_obs
    if n < 3 or n > max_cells:
        return ds
    dendros = dict(ds.adata.uns.get("dendrograms", {}))
    if force or "total" not in dendros:
        layer = ds.default_layer
        dendros["total"] = ward_manhattan(ds.layer(layer)) | {
            "layers": [layer],
            "method": "ward.D2, Manhattan",
        }
    sem = ds.semantics
    if (
        (force or "allele" not in dendros)
        and sem.get("major_cn") is Semantics.ALLELE_CN
        and sem.get("minor_cn") is Semantics.ALLELE_CN
    ):
        x = np.hstack([ds.layer("major_cn"), ds.layer("minor_cn")])
        dendros["allele"] = ward_manhattan(x) | {
            "layers": ["major_cn", "minor_cn"],
            "method": "ward.D2, Manhattan",
        }
    ds.adata.uns["dendrograms"] = dendros
    ds.validate()
    return ds
