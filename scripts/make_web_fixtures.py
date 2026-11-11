"""Fixtures for the web tests: what the Python side makes of a file, for the browser readers
(web/src/local) to be checked against.

    python scripts/make_web_fixtures.py
"""

import gzip
import json
from pathlib import Path

import numpy as np

from cnvscope.prepare import prepare
from cnvscope.readers import load
from cnvscope.readers.demo import make_demo
from cnvscope.server import build_meta

FIX = Path(__file__).parents[1] / "tests" / "fixtures"


def dump(ds, name: str, path: Path) -> None:
    meta = build_meta(ds, name)
    # layer checksums, so the browser's copy of the matrices can be checked too
    sums = {n: float(np.nansum(ds.layer(n), dtype=np.float64)) for n in ds.layer_names}
    meta["_layer_sums"] = sums
    path.write_bytes(gzip.compress(json.dumps(meta).encode(), mtime=0))


def main() -> None:
    (FIX / "h5ad").mkdir(exist_ok=True)
    ds = make_demo(n_cells=40, seed=3)
    ds.write(FIX / "h5ad" / "demo_small.h5ad")
    dump(ds, "demo_small", FIX / "h5ad" / "demo_small.meta.json.gz")
    # the ASCAT.sc export as `cnvscope serve` shows it: reader, then QC (it has a dendrogram)
    ds = prepare(load(FIX / "ascatsc_json" / "demo_export.json.gz"), cluster=False)
    dump(ds, "demo_export", FIX / "ascatsc_json" / "demo_export.meta.json.gz")


if __name__ == "__main__":
    main()
