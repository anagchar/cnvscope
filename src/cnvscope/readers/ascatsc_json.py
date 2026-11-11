"""ASCAT.sc results exported for the web: ``ascat_data.json`` (or ``.json.gz``).

This is the format of the ASCAT.sc Explorer (``R/ascatsc_to_web.R``: ``rds_to_web()``) and the
ASCAT.sc-py viewer. Fields, all keyed by cell where per cell:

- ``bins`` {chr, start, end, start_cum, end_cum}; ``chr_info``
- ``profiles`` integer total copy number; ``nMajor`` / ``nMinor`` allele copy number
- ``raw`` per-bin values: copy-number scale (R export, ploidy * 2^logR) or logR (Python
  export, has negative values)
- ``reads`` per-bin read counts; ``ci`` {nMajor_lower, nMajor_upper, nMinor_lower,
  nMinor_upper} BAF 5–95% bounds
- ``quality`` / ``qc_metrics`` per-cell QC; ``cell_types``
- ``clustering_order``, ``dendrogram`` {merge, height} (R hclust; leaves index ``profiles``)
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ..analysis import hclust_order
from ..model import CNVDataset, Modality, Semantics
from . import Reader, register


def _open(path: Path) -> dict:
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw)


def _matrix(d: dict | None, cells: list[str], nb: int) -> np.ndarray | None:
    if not d:
        return None
    m = np.full((len(cells), nb), np.nan, dtype=np.float32)
    for i, c in enumerate(cells):
        v = d.get(c)
        if v is not None:
            m[i] = np.array([np.nan if x is None else x for x in v], dtype=np.float32)
    return m


def from_dict(d: dict, name: str = "ascatsc") -> CNVDataset:
    if "bins" not in d or "profiles" not in d:
        raise ValueError("not an ASCAT.sc web export: missing 'bins' or 'profiles'")
    cells = list(d["profiles"].keys())
    b = d["bins"]
    bins = pd.DataFrame({k: b[k] for k in ("chr", "start", "end") if k in b})
    if "start_cum" in b and "end_cum" in b:
        bins["start_cum"] = np.asarray(b["start_cum"], dtype=np.int64)
        bins["end_cum"] = np.asarray(b["end_cum"], dtype=np.int64)
    nb = len(bins)

    layers: dict[str, np.ndarray] = {"total_cn": _matrix(d["profiles"], cells, nb)}
    sem: dict[str, Semantics] = {"total_cn": Semantics.INTEGER_CN}
    if d.get("nMajor") and d.get("nMinor"):
        layers["major_cn"] = _matrix(d["nMajor"], cells, nb)
        layers["minor_cn"] = _matrix(d["nMinor"], cells, nb)
        sem["major_cn"] = sem["minor_cn"] = Semantics.ALLELE_CN
        ci = d.get("ci") or {}
        if ci:
            for allele, base in (("nMajor", "major_cn"), ("nMinor", "minor_cn")):
                for side in ("lower", "upper"):
                    key = f"{allele}_{side}"
                    layers[f"{base}_{side}"] = _matrix(
                        {c: v.get(key) for c, v in ci.items() if v.get(key) is not None},
                        cells,
                        nb,
                    )
                    sem[f"{base}_{side}"] = Semantics.ALLELE_CN
    raw = _matrix(d.get("raw"), cells, nb)
    if raw is not None:
        if np.nanmin(raw[: min(20, len(cells))]) < 0:
            layers["logr"], sem["logr"] = raw, Semantics.CONTINUOUS_RATIO
        else:
            layers["cn_scale"], sem["cn_scale"] = raw, Semantics.CONTINUOUS_CN
    reads = _matrix(d.get("reads"), cells, nb)
    if reads is not None:
        layers["read_counts"], sem["read_counts"] = reads, Semantics.COUNT

    obs = pd.DataFrame(index=cells)
    for src in ("quality", "qc_metrics"):
        for c, q in (d.get(src) or {}).items():
            if c not in obs.index or not isinstance(q, dict):
                continue
            for k, v in q.items():
                if isinstance(v, int | float) and not isinstance(v, bool):
                    obs.loc[c, k] = float(v)
    ct = d.get("cell_types")
    if ct:
        obs["cell_type"] = pd.Categorical([ct.get(c) for c in cells])

    extra: dict = {}
    dend = d.get("dendrogram")
    if dend and dend.get("merge") and len(dend["merge"]) == len(cells) - 1:
        merge = np.asarray(dend["merge"], dtype=np.int64)
        order = hclust_order(merge)
        co = d.get("clustering_order")
        if co and len(co) == len(cells):
            idx = {c: i for i, c in enumerate(cells)}
            if all(isinstance(x, str) for x in co) and set(co) == set(cells):
                order = np.array([idx[c] for c in co], dtype=np.int64)
            elif all(isinstance(x, int | float) for x in co):
                o = np.asarray(co, dtype=np.int64)
                o = o - 1 if o.min() == 1 else o
                if sorted(o.tolist()) == list(range(len(cells))):
                    order = o
        extra["dendrograms"] = {
            "total": {
                "merge": merge,
                "height": np.asarray(dend["height"], dtype=np.float64),
                "order": order,
                "layers": ["total_cn"],
                "method": "ward.D2, Manhattan (from the export)",
            }
        }
    return CNVDataset.build(
        cells=obs,
        bins=bins,
        layers=layers,
        semantics=sem,
        modality=Modality.SCDNA,
        tool="ASCAT.sc",
        genome="hg38",
        default_layer="total_cn",
        extra_uns=extra | {"source": name},
    )


def _detect(path: Path) -> bool:
    n = path.name.lower()
    if not path.is_file() or not (n.endswith(".json") or n.endswith(".json.gz")):
        return False
    head = path.read_bytes()[:4096]
    if head[:2] == b"\x1f\x8b":
        with gzip.open(path, "rb") as f:
            head = f.read(4096)
    return b'"bins"' in head or b'"profiles"' in head or b'"metadata"' in head


def _read(path: Path, **_opts) -> CNVDataset:
    return from_dict(_open(path), name=path.name.split(".json")[0])


register(
    Reader(
        "ascatsc_json",
        "scDNA",
        _detect,
        _read,
        "ASCAT.sc web export (ascat_data.json[.gz], from rds_to_web())",
    )
)
