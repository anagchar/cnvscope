"""HTTP API serving one dataset to the web UI.

- ``GET  /api/meta``          the dataset: cells and their annotations/QC, bins,
                              chromosomes, layers, dendrograms. ``{"dataset": null}`` when
                              nothing is open (the UI then shows its upload screen).
- ``GET  /api/layer/{name}``  a cells × bins matrix as raw little-endian float32, rows in
                              the dataset's own order. ``?start=&stop=`` selects rows so the
                              UI can fetch large datasets in chunks. Shape is in headers.
- ``GET  /api/cell/{index}``  every layer for one cell.
- ``POST /api/open?name=``    open an uploaded file (body = file bytes): .h5ad, an ASCAT.sc
                              web export (.json/.json.gz), or any format a reader detects.
- ``POST /api/demo``          open the synthetic demo dataset.
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from fastapi import FastAPI, HTTPException, Query, Request, Response
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

from . import __version__, readers
from .model import CNVDataset
from .prepare import prepare

STATIC = Path(__file__).parent / "static"
MAX_UPLOAD = 2 * 1024**3


def _jsonable_column(s: pd.Series) -> list:
    if (
        isinstance(s.dtype, pd.CategoricalDtype)
        or s.dtype == object
        or pd.api.types.is_string_dtype(s)
    ):
        return [None if pd.isna(v) else str(v) for v in s]
    if s.dtype == bool:
        return [bool(v) for v in s]
    vals = s.to_numpy(dtype=np.float64, na_value=np.nan)
    return [None if not math.isfinite(v) else float(v) for v in vals]


def _kind(s: pd.Series) -> str:
    if s.dtype == bool:
        return "boolean"
    if pd.api.types.is_numeric_dtype(s):
        return "numeric"
    return "categorical"


def build_meta(ds: CNVDataset, name: str = "") -> dict:
    a = ds.adata
    var = a.var
    sem = ds.semantics
    return {
        "version": __version__,
        "name": name or str(a.uns.get("source", ds.tool)),
        "tool": ds.tool,
        "modality": ds.modality.value,
        "genome": a.uns.get("genome"),
        "n_cells": int(a.n_obs),
        "n_bins": int(a.n_vars),
        "cells": a.obs.index.astype(str).tolist(),
        "obs": {c: _jsonable_column(a.obs[c]) for c in a.obs.columns},
        "obs_kinds": {c: _kind(a.obs[c]) for c in a.obs.columns},
        "bins": {
            "chr": var["chr"].astype(str).tolist(),
            "start": var["start"].astype(np.int64).tolist(),
            "end": var["end"].astype(np.int64).tolist(),
            "start_cum": var["start_cum"].astype(np.int64).tolist(),
            "end_cum": var["end_cum"].astype(np.int64).tolist(),
        },
        "chr_info": ds.chr_info(),
        "layers": [{"name": n, "semantics": sem[n].value} for n in ds.track_layers],
        "intervals": {k: list(v) for k, v in ds.interval_layers.items()},
        "default_layer": ds.default_layer,
        "dendrograms": {
            k: {
                "merge": np.asarray(d["merge"]).tolist(),
                "height": np.asarray(d["height"]).tolist(),
                "order": np.asarray(d["order"]).tolist(),
                "layers": [str(x) for x in d.get("layers", [])],
                "method": str(d.get("method", "")),
            }
            for k, d in ds.dendrograms.items()
        },
    }


class State:
    def __init__(self, ds: CNVDataset | None, name: str = ""):
        self.ds: CNVDataset | None = None
        self.meta: dict | None = None
        if ds is not None:
            self.set(ds, name)

    def set(self, ds: CNVDataset, name: str = "") -> None:
        self.ds = ds
        self.meta = build_meta(ds, name)


def create_app(ds: CNVDataset | None = None, name: str = "") -> FastAPI:
    app = FastAPI(title="cnvscope", version=__version__)
    state = State(ds, name)

    def current() -> CNVDataset:
        if state.ds is None:
            raise HTTPException(404, "no dataset is open")
        return state.ds

    @app.get("/api/meta")
    def get_meta() -> dict:
        return {"dataset": state.meta}

    @app.get("/api/layer/{name}")
    def get_layer(name: str, start: int = Query(0, ge=0), stop: int | None = Query(None, ge=0)):
        ds = current()
        if name not in ds.layer_names:
            raise HTTPException(404, f"no layer {name!r}")
        n = ds.adata.n_obs
        stop = n if stop is None else min(stop, n)
        chunk = np.ascontiguousarray(ds.adata.layers[name][start:stop], dtype="<f4")
        return Response(
            chunk.tobytes(),
            media_type="application/octet-stream",
            headers={"X-Rows": str(chunk.shape[0]), "X-Cols": str(ds.adata.n_vars)},
        )

    @app.get("/api/cell/{index}")
    def get_cell(index: int) -> dict:
        ds = current()
        if not 0 <= index < ds.adata.n_obs:
            raise HTTPException(404, "no such cell")
        return {
            "cell": str(ds.adata.obs.index[index]),
            "layers": {
                k: [None if not math.isfinite(v) else float(v) for v in ds.layer(k)[index]]
                for k in ds.layer_names
            },
        }

    @app.post("/api/open")
    async def open_file(request: Request, name: str = Query(..., min_length=1)) -> dict:
        body = await request.body()
        if len(body) > MAX_UPLOAD:
            raise HTTPException(413, "file too large")
        safe = Path(name).name
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / safe
            p.write_bytes(body)
            try:
                ds = readers.load(p)
            except Exception as e:  # noqa: BLE001 - report any reader failure to the UI
                raise HTTPException(400, f"Could not open {safe}: {e}") from e
        prepare(ds)
        stem = safe.split(".json")[0].removesuffix(".h5ad").rsplit(".R", 1)[0]
        state.set(ds, stem)
        return {"dataset": state.meta}

    @app.post("/api/demo")
    def open_demo() -> dict:
        state.set(readers.load("demo"), "Demo dataset")
        return {"dataset": state.meta}

    @app.middleware("http")
    async def revalidate_page(request: Request, call_next):
        # the page must be re-checked after an update; assets have content-hashed names
        response = await call_next(request)
        if not request.url.path.startswith(("/assets/", "/api/")):
            response.headers["Cache-Control"] = "no-cache"
        return response

    if (STATIC / "index.html").exists():
        app.mount("/", StaticFiles(directory=STATIC, html=True), name="ui")
    else:

        @app.get("/", response_class=HTMLResponse)
        def no_ui() -> str:
            return (
                "<h1>cnvscope API is running</h1><p>The web UI is not built. Run "
                "<code>cd web &amp;&amp; npm install &amp;&amp; npm run build</code>, or use "
                "<code>npm run dev</code> and open http://localhost:5173.</p>"
            )

    return app
