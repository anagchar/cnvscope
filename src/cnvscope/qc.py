"""Per-cell QC metrics, as ASCAT.sc's getQCs() and the ASCAT.sc Explorer compute them.

Metrics are added as numeric ``obs`` columns. Values a reader already put there (for example
the tool's own QC table) are kept; only missing ones are computed. What can be computed
depends on the layers present:

- read counts (``count`` layer): ``n_reads_total``, ``reads_per_bin``, ``gini``
- integer copy number: ``ploidy``, ``breakpoints``, ``state_mode``, ``bins_with_cna``
- plus a log ratio or copy-number-scale track: ``mapd``, ``median_residual``, ``spikiness``,
  ``MBRSM_dispersion``, ``autocorr``, ``distance_integer_logR``, ``MSRSI_non_integerness``

Segments are runs of one integer state within a chromosome; the solution is taken at purity 1
with the profile's mean copy number as ploidy (single cells).

Everything is computed for all cells at once (numpy over the cells × bins matrices), which is
what keeps loading fast; tests/python/qc_reference.py holds the one-cell-at-a-time version
the results are checked against. The result is recorded in ``uns["qc"]`` so a saved dataset
is not recomputed when it is opened again.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from .model import CNVDataset, Semantics

QC_VERSION = 1


def gini(values: np.ndarray) -> float:
    """Gini coefficient of one cell's read counts (0 = even coverage)."""
    return float(_gini_rows(np.asarray(values, dtype=np.float64)[None, :])[0])


def _gini_rows(R: np.ndarray) -> np.ndarray:
    V = np.where(np.isfinite(R) & (R >= 0), R, np.nan)
    S = np.sort(V, axis=1)  # NaNs last, so the ranks 1..n line up with the values
    n = np.isfinite(S).sum(axis=1)
    idx = np.arange(1, S.shape[1] + 1)
    weighted = np.nansum(S * idx, axis=1)
    total = np.nansum(S, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        g = 2 * weighted / (n * total) - (n + 1) / n
    return np.where(n == 0, np.nan, np.where(total == 0, 0.0, g))


def _nanmedian_rows(X: np.ndarray) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(X, axis=1)


def _mad_rows(X: np.ndarray) -> np.ndarray:
    med = _nanmedian_rows(X)
    return 1.4826 * _nanmedian_rows(np.abs(X - med[:, None]))


def _layer_with(ds: CNVDataset, sem: Semantics, prefer: tuple[str, ...] = ()) -> str | None:
    s = ds.semantics
    for name in prefer:
        if s.get(name) is sem:
            return name
    return next((n for n in ds.track_layers if s[n] is sem), None)


def _mode_rows(P: np.ndarray) -> np.ndarray:
    """Most frequent value per row; the smallest wins ties; NaN for rows with no value."""
    fin = np.isfinite(P)
    out = np.full(P.shape[0], np.nan)
    vals = P[fin]
    if vals.size and np.all(vals == np.round(vals)) and vals.min() >= 0 and vals.max() < 10_000:
        counts = np.zeros((P.shape[0], int(vals.max()) + 1), dtype=np.int64)
        r, _ = np.nonzero(fin)
        np.add.at(counts, (r, vals.astype(np.int64)), 1)
        has = fin.any(axis=1)
        out[has] = counts[has].argmax(axis=1)
        return out
    for i in range(P.shape[0]):  # non-integer states: one cell at a time
        v = P[i][fin[i]]
        if v.size:
            u, c = np.unique(v, return_counts=True)
            out[i] = u[np.argmax(c)]
    return out


def qc_matrices(
    chrs: np.ndarray,
    width_mb: np.ndarray,
    P: np.ndarray | None,
    L: np.ndarray | None,
    R: np.ndarray | None,
    baseline: float = 2,
) -> dict[str, np.ndarray]:
    """All metrics for all cells: cells × bins matrices in, one array per metric out."""
    out: dict[str, np.ndarray] = {}
    nb = len(chrs)
    same = np.r_[False, chrs[1:] == chrs[:-1]]  # bin i continues the chromosome of bin i-1
    edges = np.flatnonzero(~same)
    slices = list(zip(edges, [*edges[1:], nb], strict=True))

    if R is not None:
        R = R.astype(np.float64)
        ok = np.isfinite(R)
        out["n_reads_total"] = np.where(ok, R, 0).sum(axis=1)
        n = ok.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            out["reads_per_bin"] = np.where(n > 0, out["n_reads_total"] / n, np.nan)
        out["gini"] = _gini_rows(R)
    if P is None:
        return out

    P = P.astype(np.float64)
    fin = np.isfinite(P)
    nfin = fin.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        psi = np.where(nfin > 0, np.where(fin, P, 0).sum(axis=1) / nfin, np.nan)
    out["ploidy"] = psi
    out["bins_with_cna"] = (fin & (np.round(P) != baseline)).sum(axis=1).astype(np.float64)
    step = fin[:, 1:] & fin[:, :-1] & (P[:, 1:] != P[:, :-1]) & same[None, 1:]
    out["breakpoints"] = step.sum(axis=1).astype(np.float64)
    out["state_mode"] = _mode_rows(P)
    if L is None:
        return out

    L = L.astype(np.float64)
    ok_cell = np.isfinite(psi) & (psi > 0)

    def nan_rows(x: np.ndarray) -> np.ndarray:
        return np.where(ok_cell, x, np.nan)

    d1 = np.where(same[None, 1:], np.abs(np.diff(L, axis=1)), np.nan)
    out["mapd"] = nan_rows(_nanmedian_rows(d1))
    with np.errstate(over="ignore", invalid="ignore"):
        out["median_residual"] = nan_rows(_nanmedian_rows(np.abs(psi[:, None] * 2.0**L - P)))
    same3 = same[1:-1] & same[2:]  # bins i, i+1, i+2 on one chromosome
    d2 = np.where(same3[None, :], np.abs(L[:, 2:] - 2 * L[:, 1:-1] + L[:, :-2]), np.nan)
    out["spikiness"] = nan_rows(_nanmedian_rows(d2))

    # segments: runs of one state within a chromosome (a missing state is a run of its own)
    n = P.shape[0]
    new = np.ones_like(fin)
    new[:, 1:] = ~same[None, 1:] | ~(P[:, 1:] == P[:, :-1])
    run = np.cumsum(new, axis=1) - 1
    gid = (run + np.arange(n)[:, None] * nb).ravel()
    finL = np.isfinite(L)
    m_runs = n * nb
    cnt = np.bincount(gid, weights=finL.ravel(), minlength=m_runs)
    sumL = np.bincount(gid, weights=np.where(finL, L, 0).ravel(), minlength=m_runs)
    runfin = np.bincount(gid, weights=fin.ravel(), minlength=m_runs) > 0
    wsum = np.bincount(gid, weights=np.broadcast_to(width_mb, P.shape).ravel(), minlength=m_runs)
    valid = (cnt > 0) & runfin
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(valid, sumL / cnt, np.nan)
    seg = mean[gid].reshape(P.shape)
    run_row = np.arange(m_runs) // nb
    with np.errstate(over="ignore", invalid="ignore"):
        cn_u = psi[run_row] * 2.0**mean
        dev = np.where(valid, wsum * (cn_u - np.round(cn_u)) ** 2, 0.0)
    dsum = np.bincount(run_row, weights=dev, minlength=n)
    dw = np.bincount(run_row, weights=np.where(valid, wsum, 0.0), minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["distance_integer_logR"] = nan_rows(np.where(dw > 0, dsum / dw, np.nan))

    resid = L - seg
    out["MBRSM_dispersion"] = nan_rows(_mad_rows(resid))

    # lag-1 autocorrelation of the residuals, per chromosome, weighted by its bins
    acs, aw = np.zeros(n), np.zeros(n)
    for a, b in slices:
        if b - a < 3:
            continue
        x, y = resid[:, a : b - 1], resid[:, a + 1 : b]
        ok = np.isfinite(x) & np.isfinite(y)
        k = ok.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            mx = np.where(ok, x, 0.0).sum(1) / k
            my = np.where(ok, y, 0.0).sum(1) / k
            dx = np.where(ok, x - mx[:, None], 0.0)
            dy = np.where(ok, y - my[:, None], 0.0)
            sxx, syy, sxy = (dx * dx).sum(1), (dy * dy).sum(1), (dx * dy).sum(1)
            c = sxy / np.sqrt(sxx * syy)
        use = (k >= 2) & (sxx > 0) & (syy > 0) & np.isfinite(c)
        acs += np.where(use, c * (b - a), 0.0)
        aw += np.where(use, b - a, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out["autocorr"] = nan_rows(np.where(aw > 0, acs / aw, np.nan))

    # distance of segment logR from the logR its integer state predicts, by (chromosome, state)
    chr_idx = np.cumsum(~same) - 1
    use = fin & (P > 0) & np.isfinite(seg) & ok_cell[:, None]
    r, c = np.nonzero(use)
    out["MSRSI_non_integerness"] = np.full(n, np.nan)
    if r.size:
        with np.errstate(divide="ignore", invalid="ignore"):
            val = seg[r, c] - np.log2(P[r, c] / psi[r])
        g = pd.DataFrame({"row": r, "chr": chr_idx[c], "st": P[r, c], "v": val})
        med = g.groupby(["row", "chr", "st"], sort=False)["v"].median().rename("m").reset_index()
        med["d"] = (med["m"] - med.groupby("row")["m"].transform("median")).abs()
        mad = 1.4826 * med.groupby("row")["d"].median()
        out["MSRSI_non_integerness"][mad.index.to_numpy()] = mad.to_numpy()
    return out


def compute_qc(ds: CNVDataset, force: bool = False) -> CNVDataset:
    """Fill in QC metrics as ``obs`` columns (existing columns are kept)."""
    a = ds.adata
    done = a.uns.get("qc")
    if not force and isinstance(done, dict) and int(done.get("version", 0)) == QC_VERSION:
        computed = [str(k) for k in np.atleast_1d(done.get("computed", []))]
        if all(k in a.obs.columns for k in computed):
            return ds
    var = a.var
    chrs = var["chr"].astype(str).to_numpy()
    width_mb = (var["end"].to_numpy() - var["start"].to_numpy()) / 1e6

    total = _layer_with(ds, Semantics.INTEGER_CN, ("total_cn",))
    reads = _layer_with(ds, Semantics.COUNT, ("read_counts",))
    ratio = _layer_with(ds, Semantics.CONTINUOUS_RATIO, ("logr",))
    cn_scale = _layer_with(ds, Semantics.CONTINUOUS_CN)
    P = ds.layer(total) if total else None
    R = ds.layer(reads) if reads else None
    L = ds.layer(ratio) if ratio else None
    if L is None and cn_scale and P is not None:
        # bins on the copy-number scale → logR around the cell's mean copy number
        fin = np.isfinite(P)
        k = fin.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            psi = np.where(k > 0, np.where(fin, P, 0).sum(axis=1) / k, np.nan)
            CS = ds.layer(cn_scale)
            L = np.where(CS > 0, np.log2(CS / psi[:, None]), np.nan)
    if P is None and R is None:
        return ds

    metrics = qc_matrices(chrs, width_mb, P, L, R)
    computed = []
    for k, col in metrics.items():
        col = np.asarray(col, dtype=np.float64)
        if not np.isfinite(col).any():
            continue  # a metric no cell has is left out
        if k in a.obs.columns:
            have = a.obs[k].to_numpy(dtype=np.float64, na_value=np.nan)
            col = np.where(np.isfinite(have), have, col)
        a.obs[k] = col
        computed.append(k)
    a.uns["qc"] = {"version": QC_VERSION, "computed": sorted(computed)}
    return ds
