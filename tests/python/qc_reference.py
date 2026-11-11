"""Reference: the original per-cell QC implementation (one cell at a time).

The vectorised cnvscope.qc must give the same numbers; test_qc_fast.py compares them.
"""

# ruff: noqa
from __future__ import annotations

import numpy as np

from cnvscope.model import CNVDataset, Semantics


def gini(values: np.ndarray) -> float:
    v = np.sort(values[np.isfinite(values) & (values >= 0)])
    n = len(v)
    if n == 0:
        return np.nan
    s = v.sum()
    if s == 0:
        return 0.0
    return float(2 * np.sum(np.arange(1, n + 1) * v) / (n * s) - (n + 1) / n)


def _mad_scaled(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    if not len(x):
        return np.nan
    m = np.median(x)
    return float(1.4826 * np.median(np.abs(x - m)))


def _median(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(np.median(x)) if len(x) else np.nan


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    da, db = a - a.mean(), b - b.mean()
    saa, sbb = (da * da).sum(), (db * db).sum()
    return float((da * db).sum() / np.sqrt(saa * sbb)) if saa > 0 and sbb > 0 else np.nan


def _layer_with(ds: CNVDataset, sem: Semantics, prefer: tuple[str, ...] = ()) -> str | None:
    s = ds.semantics
    for name in prefer:
        if s.get(name) is sem:
            return name
    return next((n for n in ds.track_layers if s[n] is sem), None)


def cell_metrics(
    chr_slices: list[tuple[int, int]],
    width_mb: np.ndarray,
    chrs: np.ndarray,
    prof: np.ndarray | None,
    logr: np.ndarray | None,
    reads: np.ndarray | None,
    baseline: float = 2,
) -> dict[str, float]:
    r: dict[str, float] = {}
    if reads is not None:
        ok = np.isfinite(reads)
        r["n_reads_total"] = float(reads[ok].sum())
        r["reads_per_bin"] = float(reads[ok].mean()) if ok.any() else np.nan
        r["gini"] = gini(reads)
    if prof is None:
        return r
    nb = len(prof)
    fin = np.isfinite(prof)
    psi = float(prof[fin].mean()) if fin.any() else np.nan
    r["ploidy"] = psi
    r["bins_with_cna"] = float(np.sum(fin & (np.round(prof) != baseline)))
    bp = 0
    for a, b in chr_slices:
        p = prof[a:b]
        both = np.isfinite(p[1:]) & np.isfinite(p[:-1])
        bp += int(np.sum(both & (p[1:] != p[:-1])))
    r["breakpoints"] = float(bp)
    if fin.any():
        vals, counts = np.unique(prof[fin], return_counts=True)
        r["state_mode"] = float(vals[np.argmax(counts)])  # smallest state wins ties
    if logr is None or not np.isfinite(psi) or psi <= 0:
        return r

    L = logr.astype(np.float64)
    d1, d2 = [], []
    for a, b in chr_slices:
        x = L[a:b]
        d1.append(np.abs(np.diff(x)))
        if b - a > 2:
            d2.append(np.abs(x[2:] - 2 * x[1:-1] + x[:-2]))
    r["mapd"] = _median(np.concatenate(d1)) if d1 else np.nan
    with np.errstate(over="ignore", invalid="ignore"):
        r["median_residual"] = _median(np.abs(psi * 2.0**L - prof))
    r["spikiness"] = _median(np.concatenate(d2)) if d2 else np.nan

    seg = np.full(nb, np.nan)
    dw = dsum = 0.0
    for a, b in chr_slices:
        s = a
        while s < b:
            e = s + 1
            while e < b and prof[e] == prof[s]:
                e += 1
            x = L[s:e]
            okx = np.isfinite(x)
            if okx.any() and np.isfinite(prof[s]):
                m = x[okx].mean()
                seg[s:e] = m
                cn_u = psi * 2.0**m
                w = width_mb[s:e].sum()
                dsum += w * (cn_u - round(cn_u)) ** 2
                dw += w
            s = e
    if dw > 0:
        r["distance_integer_logR"] = dsum / dw
    resid = L - seg
    r["MBRSM_dispersion"] = _mad_scaled(resid)
    acs = aw = 0.0
    for a, b in chr_slices:
        if b - a < 3:
            continue
        x, y = resid[a : b - 1], resid[a + 1 : b]
        ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 2:
            continue
        c = _corr(x[ok], y[ok])
        if np.isfinite(c):
            acs += c * (b - a)
            aw += b - a
    if aw:
        r["autocorr"] = acs / aw
    groups: dict[tuple[str, float], list[float]] = {}
    for i in range(nb):
        st = prof[i]
        if not (np.isfinite(st) and st > 0) or not np.isfinite(seg[i]):
            continue
        groups.setdefault((chrs[i], st), []).append(seg[i] - np.log2(st / psi))
    if groups:
        r["MSRSI_non_integerness"] = _mad_scaled(np.array([np.median(v) for v in groups.values()]))
    return r


def compute_qc(ds: CNVDataset) -> CNVDataset:
    """Fill in QC metrics as ``obs`` columns (existing columns are kept)."""
    a = ds.adata
    var = a.var
    chrs = var["chr"].astype(str).to_numpy()
    edges = np.flatnonzero(chrs[1:] != chrs[:-1]) + 1
    bounds = [0, *edges.tolist(), len(chrs)]
    slices = list(zip(bounds[:-1], bounds[1:], strict=True))
    width_mb = (var["end"].to_numpy() - var["start"].to_numpy()) / 1e6

    total = _layer_with(ds, Semantics.INTEGER_CN, ("total_cn",))
    reads = _layer_with(ds, Semantics.COUNT, ("read_counts",))
    ratio = _layer_with(ds, Semantics.CONTINUOUS_RATIO, ("logr",))
    cn_scale = _layer_with(ds, Semantics.CONTINUOUS_CN)
    P = ds.layer(total) if total else None
    R = ds.layer(reads) if reads else None
    LR = ds.layer(ratio) if ratio else None
    CS = ds.layer(cn_scale) if cn_scale and LR is None else None
    if P is None and R is None:
        return ds

    rows = []
    for i in range(a.n_obs):
        prof = P[i] if P is not None else None
        logr = LR[i] if LR is not None else None
        if logr is None and CS is not None and prof is not None:
            fin = np.isfinite(prof)
            psi = prof[fin].mean() if fin.any() else np.nan
            with np.errstate(divide="ignore", invalid="ignore"):
                logr = np.where(CS[i] > 0, np.log2(CS[i] / psi), np.nan)
        rows.append(
            cell_metrics(slices, width_mb, chrs, prof, logr, R[i] if R is not None else None)
        )
    keys = sorted({k for r in rows for k in r})
    for k in keys:
        col = np.array([r.get(k, np.nan) for r in rows], dtype=np.float64)
        if k in a.obs.columns:
            have = a.obs[k].to_numpy(dtype=np.float64, na_value=np.nan)
            col = np.where(np.isfinite(have), have, col)
        a.obs[k] = col
    return ds
