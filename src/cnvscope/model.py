"""The canonical dataset every reader produces and the UI consumes.

A dataset is an AnnData object:

- ``obs``    one row per cell (index = barcode); optional columns such as ``cluster``,
             ``clone``, ``cell_type``, ``sample``, ``ploidy``, QC metrics, ``passing``.
- ``var``    one row per genomic bin: ``chr``, ``start``, ``end`` (1-based inclusive or
             0-based half-open, as the tool gives them) and ``start_cum``/``end_cum`` on the
             concatenated genome.
- ``layers`` cells × bins value tracks (float32, NaN = missing), e.g. ``total_cn``,
             ``major_cn``, ``minor_cn``, ``logr``, ``rdr``, ``baf``, ``state``.
- ``uns``    ``modality``, ``tool``, ``genome``, ``semantics`` (layer → :class:`Semantics`),
             ``default_layer`` and ``dendrograms``: name → {``merge`` (R hclust layout,
             n-1 × 2, negative = leaf), ``height``, ``order`` (leaf order, row indices),
             ``layers``, ``method``}. ``total`` clusters the default layer; ``allele``
             clusters major and minor copy number together.

Layers named ``<layer>_lower`` / ``<layer>_upper`` are an interval around ``<layer>``
(e.g. BAF-derived 5–95% bounds of the major allele) and are not shown as tracks of their own.

The UI decides colours, baseline and "is this altered" from ``semantics`` and ``modality``,
never from the tool name.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


class Modality(str, Enum):
    SCDNA = "scDNA"
    SCATAC = "scATAC"
    SCRNA = "scRNA"
    MULTIOME = "multiome"


class Semantics(str, Enum):
    INTEGER_CN = "integer_cn"  # absolute copy number 0, 1, 2, ... (baseline = ploidy, usually 2)
    ALLELE_CN = "allele_cn"  # copy number of one allele/haplotype (major, minor)
    TERNARY_STATE = "ternary_state"  # -1 loss, 0 neutral, +1 gain
    CONTINUOUS_RATIO = "continuous_ratio"  # log2 ratio, centred on 0
    CONTINUOUS_CN = "continuous_cn"  # bin values on the copy-number scale (ploidy * 2^logR)
    RATIO = "ratio"  # linear ratio, centred on 1 (RDR, expression ratio)
    BAF = "baf"  # B-allele frequency, 0..1, balanced at 0.5
    PROBABILITY = "probability"  # posterior 0..1
    COUNT = "count"  # read counts and other non-negative counts


REQUIRED_VAR = ("chr", "start", "end")

# Natural order of chromosome names; anything unknown sorts after Y/M by name.
_CHR_RANK = {str(i): i for i in range(1, 23)} | {"X": 23, "Y": 24, "M": 25, "MT": 25}


def chr_sort_key(name: str) -> tuple[int, str]:
    bare = str(name).removeprefix("chr")
    return (_CHR_RANK.get(bare, 100), bare)


class ValidationError(ValueError):
    pass


@dataclass
class CNVDataset:
    adata: ad.AnnData

    # ---------- construction ----------
    @classmethod
    def build(
        cls,
        *,
        cells: pd.DataFrame | list[str],
        bins: pd.DataFrame,
        layers: dict[str, np.ndarray],
        semantics: dict[str, Semantics | str],
        modality: Modality | str,
        tool: str,
        genome: str = "hg38",
        default_layer: str | None = None,
        extra_uns: dict | None = None,
    ) -> CNVDataset:
        obs = cells if isinstance(cells, pd.DataFrame) else pd.DataFrame(index=list(cells))
        obs.index = obs.index.astype(str)
        var = bins.reset_index(drop=True).copy()
        var.index = [f"bin{i}" for i in range(len(var))]
        if not {"start_cum", "end_cum"} <= set(var.columns):
            var = add_cumulative(var)
        first = next(iter(layers))
        adata = ad.AnnData(
            X=None,  # every track is a named layer; X would only duplicate one
            obs=obs,
            var=var,
            layers={k: np.asarray(v, dtype=np.float32) for k, v in layers.items()},
        )
        adata.uns.update(
            modality=Modality(modality).value,
            tool=tool,
            genome=genome,
            semantics={k: Semantics(v).value for k, v in semantics.items()},
            default_layer=default_layer or first,
        )
        if extra_uns:
            adata.uns.update(extra_uns)
        ds = cls(adata)
        ds.validate()
        return ds

    # ---------- accessors ----------
    @property
    def modality(self) -> Modality:
        return Modality(self.adata.uns["modality"])

    @property
    def tool(self) -> str:
        return str(self.adata.uns["tool"])

    @property
    def semantics(self) -> dict[str, Semantics]:
        return {k: Semantics(v) for k, v in self.adata.uns["semantics"].items()}

    @property
    def default_layer(self) -> str:
        return str(self.adata.uns["default_layer"])

    @property
    def dendrograms(self) -> dict[str, dict]:
        return dict(self.adata.uns.get("dendrograms", {}))

    @property
    def order(self) -> np.ndarray:
        """Row indices of cells in display order (the total-CN dendrogram's leaves if any)."""
        d = self.dendrograms
        first = d.get("total") or next(iter(d.values()), None)
        if first is None:
            return np.arange(self.adata.n_obs)
        return np.asarray(first["order"], dtype=np.int64)

    @property
    def interval_layers(self) -> dict[str, tuple[str, str]]:
        """Layer → (lower, upper) layers bounding it."""
        names = set(self.layer_names)
        return {
            n: (f"{n}_lower", f"{n}_upper")
            for n in names
            if f"{n}_lower" in names and f"{n}_upper" in names
        }

    @property
    def track_layers(self) -> list[str]:
        """Layers shown as tracks: all but interval bounds."""
        bounds = {b for pair in self.interval_layers.values() for b in pair}
        return [n for n in self.layer_names if n not in bounds]

    @property
    def layer_names(self) -> list[str]:
        # anndata >= 0.12 lists X as layers[None]; only named layers are tracks
        return [k for k in self.adata.layers.keys() if k is not None]

    def layer(self, name: str) -> np.ndarray:
        return np.asarray(self.adata.layers[name], dtype=np.float32)

    def chr_info(self) -> list[dict]:
        var = self.adata.var
        out = []
        for c, idx in var.groupby("chr", sort=False).indices.items():
            a, b = int(idx.min()), int(idx.max()) + 1
            out.append(
                {
                    "chr": str(c),
                    "bin_start": a,
                    "bin_end": b,
                    "start_cum": int(var["start_cum"].iloc[a]),
                    "end_cum": int(var["end_cum"].iloc[b - 1]),
                }
            )
        return out

    # ---------- validation ----------
    def validate(self) -> None:
        a = self.adata
        for col in REQUIRED_VAR:
            if col not in a.var:
                raise ValidationError(f"bins are missing column {col!r}")
        if not a.obs.index.is_unique:
            raise ValidationError("cell barcodes are not unique")
        sem = a.uns.get("semantics", {})
        for name in self.layer_names:
            if name not in sem:
                raise ValidationError(f"layer {name!r} has no semantics")
            Semantics(sem[name])
        for name in sem:
            if name not in self.layer_names:
                raise ValidationError(f"semantics given for missing layer {name!r}")
        if a.uns.get("default_layer") not in self.layer_names:
            raise ValidationError("default_layer is not a layer")
        Modality(a.uns.get("modality"))
        # bins must be grouped by chromosome and sorted within it
        chrs = a.var["chr"].astype(str).to_numpy()
        seen: set[str] = set()
        for i, c in enumerate(chrs):
            if i and c != chrs[i - 1]:
                if c in seen:
                    raise ValidationError(f"bins of {c} are not contiguous")
            seen.add(c)
        starts = a.var["start"].to_numpy()
        same = chrs[1:] == chrs[:-1]
        if np.any(np.diff(starts)[same] < 0):
            raise ValidationError("bins are not sorted by start within a chromosome")
        for name, d in a.uns.get("dendrograms", {}).items():
            o = np.asarray(d["order"])
            if sorted(o.tolist()) != list(range(a.n_obs)):
                raise ValidationError(
                    f"dendrogram {name!r}: order is not a permutation of the cells"
                )
            if np.asarray(d["merge"]).shape != (a.n_obs - 1, 2):
                raise ValidationError(f"dendrogram {name!r}: merge must have n_cells - 1 rows")

    # ---------- I/O ----------
    def write(self, path: str | Path) -> None:
        self.adata.write_h5ad(path, compression="gzip")

    @classmethod
    def read(cls, path: str | Path) -> CNVDataset:
        ds = cls(ad.read_h5ad(path))
        ds.validate()
        return ds


def add_cumulative(var: pd.DataFrame) -> pd.DataFrame:
    """Add ``start_cum``/``end_cum``: bins on the concatenated genome, chromosomes back to back.

    Each chromosome spans from 0 to its last bin's end, so gaps before the first bin are kept.
    """
    var = var.copy()
    var["chr"] = var["chr"].astype(str)
    offset = 0
    start_cum = np.zeros(len(var), dtype=np.int64)
    end_cum = np.zeros(len(var), dtype=np.int64)
    for _, idx in var.groupby("chr", sort=False).indices.items():
        s = var["start"].to_numpy()[idx].astype(np.int64)
        e = var["end"].to_numpy()[idx].astype(np.int64)
        start_cum[idx] = s + offset
        end_cum[idx] = e + offset
        offset += int(e.max())
    var["start_cum"] = start_cum
    var["end_cum"] = end_cum
    return var
