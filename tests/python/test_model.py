import numpy as np
import pandas as pd
import pytest

from cnvscope.analysis import scipy_to_hclust
from cnvscope.genome import segments_to_bins, uniform_bins
from cnvscope.model import CNVDataset, Semantics, ValidationError
from cnvscope.readers import REGISTRY, load
from cnvscope.readers.demo import make_demo


def tiny(**kw) -> CNVDataset:
    bins = pd.DataFrame(
        {
            "chr": ["chr1"] * 3 + ["chr2"] * 2,
            "start": [0, 10, 20, 0, 10],
            "end": [10, 20, 30, 10, 20],
        }
    )
    cn = np.array([[2, 2, 3, 2, 1], [2, 2, 2, 2, 2]], dtype=float)
    args = dict(
        cells=["a", "b"],
        bins=bins,
        layers={"total_cn": cn},
        semantics={"total_cn": "integer_cn"},
        modality="scDNA",
        tool="t",
    )
    args.update(kw)
    return CNVDataset.build(**args)


def test_build_adds_cumulative_coordinates():
    ds = tiny()
    assert ds.adata.var["start_cum"].tolist() == [0, 10, 20, 30, 40]
    assert ds.adata.var["end_cum"].tolist() == [10, 20, 30, 40, 50]
    info = ds.chr_info()
    assert [(c["chr"], c["bin_start"], c["bin_end"]) for c in info] == [
        ("chr1", 0, 3),
        ("chr2", 3, 5),
    ]


def test_semantics_are_required_for_every_layer():
    with pytest.raises(ValidationError):
        tiny(semantics={})
    with pytest.raises(ValueError):
        tiny(semantics={"total_cn": "not_a_semantic"})


def test_unsorted_bins_rejected():
    bins = pd.DataFrame({"chr": ["chr1", "chr2", "chr1"], "start": [0, 0, 10], "end": [10, 10, 20]})
    with pytest.raises(ValidationError):
        tiny(bins=bins, layers={"total_cn": np.zeros((2, 3))})


def test_h5ad_round_trip(tmp_path):
    ds = make_demo(n_cells=40)
    p = tmp_path / "d.h5ad"
    ds.write(p)
    back = load(p)
    assert back.semantics == ds.semantics
    assert back.semantics["logr"] is Semantics.CONTINUOUS_RATIO
    assert back.modality == ds.modality
    np.testing.assert_array_equal(back.order, ds.order)
    np.testing.assert_allclose(back.layer("logr"), ds.layer("logr"), equal_nan=True)
    assert back.dendrograms["total"]["merge"].shape == (39, 2)
    assert set(back.dendrograms) == {"total", "allele"}
    assert back.interval_layers["major_cn"] == ("major_cn_lower", "major_cn_upper")
    assert "major_cn_lower" not in back.track_layers


def test_scipy_to_hclust_leaves_and_merges():
    # three leaves: (0,1) at h=1, then (2, merge0) at h=2
    Z = np.array([[0, 1, 1.0, 2], [2, 3, 2.0, 3]])
    h = scipy_to_hclust(Z)
    assert h["merge"].tolist() == [[-1, -2], [-3, 1]]
    assert h["height"].tolist() == [1.0, 2.0]


def test_segments_to_bins_takes_segment_at_bin_midpoint():
    bins = pd.DataFrame({"chr": ["chr1"] * 4, "start": [0, 10, 20, 30], "end": [10, 20, 30, 40]})
    seg = pd.DataFrame(
        {
            "cell": ["x", "x", "y"],
            "chr": ["1", "1", "chr1"],
            "start": [0, 18, 0],
            "end": [17, 40, 12],
            "cn": [2, 4, 1],
        }
    )
    out = segments_to_bins(seg, bins, ["cn"], cells=["x", "y"])["cn"]
    np.testing.assert_array_equal(out[0], [2, 2, 4, 4])
    np.testing.assert_array_equal(out[1], [1, np.nan, np.nan, np.nan])


def test_uniform_bins_cover_genome():
    b = uniform_bins("hg38", 50_000_000)
    assert b["chr"].iloc[0] == "chr1" and b["chr"].iloc[-1] == "chrX"
    assert b[b.chr == "chr1"]["end"].max() == 248956422


def test_demo_registered():
    assert "demo" in REGISTRY
    ds = load("demo")
    assert ds.adata.n_obs == 300 and "total" in ds.dendrograms
    for k in ("ploidy", "gini", "mapd", "median_residual", "bins_with_cna", "n_reads_total"):
        assert k in ds.adata.obs and ds.adata.obs[k].notna().all(), k


def test_ward_d2_on_manhattan_matches_r():
    # ward.D2 Lance-Williams by hand: {0,1} at 1, then +2 at sqrt(25/3), then +3 at sqrt(98/3)
    x = np.array([[0, 0], [1, 0], [3, 0], [6, 0.0]])
    from cnvscope.analysis import ward_manhattan

    h = ward_manhattan(x)
    np.testing.assert_allclose(sorted(h["height"]), [1, np.sqrt(25 / 3), np.sqrt(98 / 3)])


def test_hclust_order_left_first():
    from cnvscope.analysis import hclust_order

    merge = np.array([[-2, -3], [-1, 1]])
    assert hclust_order(merge).tolist() == [0, 1, 2]
