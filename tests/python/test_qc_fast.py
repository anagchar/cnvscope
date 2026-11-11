"""The vectorised QC metrics give the same numbers as the one-cell-at-a-time reference."""

from pathlib import Path

import numpy as np
import pytest
import qc_reference as ref

from cnvscope.qc import qc_matrices
from cnvscope.readers import load
from cnvscope.readers.demo import make_demo

FIX = Path(__file__).parents[1] / "fixtures" / "ascatsc_json" / "demo_export.json.gz"


def reference(chrs, width_mb, P, L, R):
    edges = np.flatnonzero(chrs[1:] != chrs[:-1]) + 1
    bounds = [0, *edges.tolist(), len(chrs)]
    slices = list(zip(bounds[:-1], bounds[1:], strict=True))
    rows = [
        ref.cell_metrics(
            slices,
            width_mb,
            chrs,
            P[i] if P is not None else None,
            L[i] if L is not None else None,
            R[i] if R is not None else None,
        )
        for i in range(P.shape[0] if P is not None else R.shape[0])
    ]
    return {k: np.array([r.get(k, np.nan) for r in rows]) for k in {k for r in rows for k in r}}


def compare(chrs, width_mb, P, L, R):
    # both in float64, so the algorithms are compared rather than float32 rounding
    P, L, R = (None if M is None else np.asarray(M, dtype=np.float64) for M in (P, L, R))
    want = reference(chrs, width_mb, P, L, R)
    got = qc_matrices(chrs, width_mb, P, L, R)
    for k, v in want.items():
        np.testing.assert_allclose(got[k], v, rtol=1e-9, atol=1e-12, equal_nan=True, err_msg=k)


@pytest.mark.parametrize("seed", [0, 1])
def test_matches_reference_on_demo(seed):
    ds = make_demo(n_cells=60, seed=seed)
    v = ds.adata.var
    chrs = v["chr"].astype(str).to_numpy()
    width = (v["end"].to_numpy() - v["start"].to_numpy()) / 1e6
    P, L, R = ds.layer("total_cn"), ds.layer("logr"), ds.layer("read_counts")
    # missing values in every track, a cell with no calls and a cell with no reads
    rng = np.random.default_rng(seed)
    for M in (P, L, R):
        M[rng.random(M.shape) < 0.01] = np.nan
    P[3] = np.nan
    R[4] = np.nan
    compare(chrs, width, P, L, R)


def test_matches_reference_on_ascatsc_export():
    ds = load(FIX)
    v = ds.adata.var
    chrs = v["chr"].astype(str).to_numpy()
    width = (v["end"].to_numpy() - v["start"].to_numpy()) / 1e6
    P = ds.layer("total_cn")
    with np.errstate(divide="ignore", invalid="ignore"):
        psi = np.nanmean(P, axis=1, keepdims=True)
        L = np.where(ds.layer("cn_scale") > 0, np.log2(ds.layer("cn_scale") / psi), np.nan)
    compare(chrs, width, P, L, ds.layer("read_counts"))


def test_saved_qc_is_not_recomputed(tmp_path):
    from cnvscope.prepare import prepare

    ds = make_demo(n_cells=20)
    assert ds.adata.uns["qc"]["version"] == 1
    # blanked after QC: a recompute would fill it in again, a skip leaves it blank
    ds.adata.obs["mapd"] = np.nan
    p = tmp_path / "d.h5ad"
    ds.write(p)
    back = prepare(load(p))
    assert back.adata.obs["mapd"].isna().all()
    assert np.isfinite(back.adata.obs["gini"]).all()
