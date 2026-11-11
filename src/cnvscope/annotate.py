"""Join a table of cell annotations (cell type, sample, clone, ...) onto a dataset's cells."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .model import CNVDataset

#: column names taken as the cell barcode, best first (else the first column)
CELL_COLS = ("cell", "barcode", "cell_id", "cellid", "cell_barcode", "barcodes")
#: renamed so the viewer finds them: epiAneufinder's annot_dt uses `annot`
RENAME = {"annot": "cell_type", "celltype": "cell_type", "cell.type": "cell_type"}


def _key(name: str) -> str:
    """A barcode without the bits tools add: `cell-` prefixes, `-1` suffixes, sample prefixes."""
    s = re.sub(r"^cell[-_.]", "", str(name))
    s = re.sub(r"[-_.]\d+$", "", s)
    return s.rsplit("#", 1)[-1].rsplit(":", 1)[-1]  # ArchR "sample#BARCODE", "sample:BARCODE"


def read_table(path: Path) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    t = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False, na_values=["", "NA"])
    low = {c.lower(): c for c in t.columns}
    cell = next((low[c] for c in CELL_COLS if c in low), t.columns[0])
    t = t.set_index(cell)
    t.index = t.index.astype(str)
    return t.rename(columns={c: RENAME.get(c.lower(), c) for c in t.columns})


def annotate(ds: CNVDataset, path: str | Path) -> int:
    """Add the table's columns to ``ds``'s cells; returns how many cells were matched.

    Barcodes are matched exactly, else after removing ``cell-`` prefixes, ``-1`` suffixes and
    ``sample#`` prefixes on both sides.
    """
    t = read_table(Path(path))
    obs = ds.adata.obs
    idx = obs.index.astype(str)
    hit = idx.isin(t.index)
    if hit.sum() < len(idx):
        keys = pd.Series([_key(i) for i in t.index], index=t.index)
        if keys.duplicated().any():
            keys = keys[~keys.duplicated(keep=False)]
        by_key = pd.Series(keys.index, index=keys.to_numpy())
        mapped = pd.Index([i if i in t.index else by_key.get(_key(i)) for i in idx])
    else:
        mapped = pd.Index(idx)
    matched = int(pd.Series(mapped).notna().sum())
    if matched == 0:
        raise ValueError(
            f"no barcode in {Path(path).name} matches the dataset's cells "
            f"(e.g. {idx[0]!r} vs {t.index[0]!r})"
        )
    for c in t.columns:
        vals = t[c].reindex(mapped).to_numpy()
        num = pd.to_numeric(pd.Series(vals), errors="coerce")
        col = num.to_numpy() if num.notna().sum() == pd.Series(vals).notna().sum() else vals
        obs[c] = pd.Categorical(col) if col.dtype == object else col
    return matched
