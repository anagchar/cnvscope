from pathlib import Path

import numpy as np
from fastapi.testclient import TestClient

from cnvscope.readers.demo import make_demo
from cnvscope.server import create_app

FIX = Path(__file__).parents[1] / "fixtures" / "ascatsc_json" / "demo_export.json.gz"


def test_meta_and_layer_chunks():
    ds = make_demo(n_cells=30)
    client = TestClient(create_app(ds))
    meta = client.get("/api/meta").json()["dataset"]
    assert meta["n_cells"] == 30 and meta["modality"] == "scDNA"
    names = {x["name"] for x in meta["layers"]}
    assert {"total_cn", "major_cn", "minor_cn", "logr", "read_counts"} == names
    assert meta["intervals"]["minor_cn"] == ["minor_cn_lower", "minor_cn_upper"]
    assert len(meta["dendrograms"]["total"]["merge"]) == 29
    assert len(meta["bins"]["start_cum"]) == meta["n_bins"]
    assert meta["obs_kinds"]["cell_type"] == "categorical"
    assert meta["obs_kinds"]["gini"] == "numeric"

    r = client.get("/api/layer/total_cn", params={"start": 5, "stop": 12})
    assert r.status_code == 200
    rows, cols = int(r.headers["x-rows"]), int(r.headers["x-cols"])
    assert (rows, cols) == (7, meta["n_bins"])
    m = np.frombuffer(r.content, dtype="<f4").reshape(rows, cols)
    np.testing.assert_array_equal(m, ds.layer("total_cn")[5:12])

    assert client.get("/api/layer/nope").status_code == 404
    cell = client.get("/api/cell/0").json()
    assert cell["cell"] == meta["cells"][0] and len(cell["layers"]["logr"]) == meta["n_bins"]


def test_empty_then_open_upload_and_demo():
    client = TestClient(create_app())
    assert client.get("/api/meta").json() == {"dataset": None}
    assert client.get("/api/layer/total_cn").status_code == 404

    r = client.post("/api/open", params={"name": FIX.name}, content=FIX.read_bytes())
    assert r.status_code == 200, r.text
    meta = r.json()["dataset"]
    assert meta["tool"] == "ASCAT.sc" and meta["name"] == "demo_export"
    assert set(meta["dendrograms"]) == {"total", "allele"}
    assert client.get("/api/meta").json()["dataset"]["n_cells"] == meta["n_cells"]

    bad = client.post("/api/open", params={"name": "x.rds"}, content=b"xx")
    assert bad.status_code == 400 and "x.rds" in bad.json()["detail"]
    junk = client.post("/api/open", params={"name": "x.txt"}, content=b"hello")
    assert junk.status_code == 400

    assert client.post("/api/demo").json()["dataset"]["tool"] == "demo"


def test_page_is_revalidated():
    client = TestClient(create_app())
    assert client.get("/").headers.get("cache-control") == "no-cache"
