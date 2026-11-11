"""The ASCAT.sc web export (fixture written by the ASCAT.sc Explorer's demo generator)."""

import gzip
import json
from pathlib import Path

import numpy as np

from cnvscope.model import Semantics
from cnvscope.prepare import prepare
from cnvscope.readers import detect, load

FIX = Path(__file__).parents[1] / "fixtures" / "ascatsc_json" / "demo_export.json.gz"


def raw():
    return json.loads(gzip.decompress(FIX.read_bytes()))


def test_detected_and_read():
    assert detect(FIX).name == "ascatsc_json"
    ds = load(FIX)
    d = raw()
    cells = list(d["profiles"])
    assert ds.adata.obs.index.tolist() == cells
    assert ds.adata.n_vars == len(d["bins"]["chr"])
    sem = ds.semantics
    assert sem["total_cn"] is Semantics.INTEGER_CN
    assert sem["major_cn"] is Semantics.ALLELE_CN
    assert sem["read_counts"] is Semantics.COUNT
    # raw on the copy-number scale (no negative values) → continuous_cn
    assert sem["cn_scale"] is Semantics.CONTINUOUS_CN
    c = cells[3]
    np.testing.assert_array_equal(ds.layer("total_cn")[3], d["profiles"][c])
    np.testing.assert_array_equal(ds.layer("minor_cn")[3], d["nMinor"][c])
    # tool's own coordinates kept
    assert ds.adata.var["start_cum"].tolist() == d["bins"]["start_cum"]


def test_intervals_quality_and_cell_types():
    ds = load(FIX)
    d = raw()
    c0 = list(d["profiles"])[0]
    assert ds.interval_layers["major_cn"] == ("major_cn_lower", "major_cn_upper")
    np.testing.assert_allclose(
        ds.layer("major_cn_lower")[0], d["ci"][c0]["nMajor_lower"], rtol=1e-6
    )
    assert np.isnan(ds.layer("major_cn_lower")[5]).all()  # no CI for that cell
    assert ds.adata.obs.loc[c0, "median_residual"] == d["quality"][c0]["median_residual"]
    assert ds.adata.obs.loc[c0, "cell_type"] == d["cell_types"][c0]


def test_exported_dendrogram_kept_and_allele_one_added():
    ds = prepare(load(FIX))
    d = raw()
    cells = list(d["profiles"])
    t = ds.dendrograms["total"]
    assert np.asarray(t["merge"]).tolist() == d["dendrogram"]["merge"]
    assert [cells[i] for i in t["order"]] == d["clustering_order"]
    assert "allele" in ds.dendrograms
    # exported QC values win over computed ones; missing ones are filled in
    c0 = cells[0]
    assert ds.adata.obs.loc[c0, "mapd"] == d["quality"][c0]["mapd"]
    assert np.isfinite(ds.adata.obs["spikiness"]).all()
