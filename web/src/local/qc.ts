// Per-cell QC metrics in the browser, for files opened without the cnvscope
// server. A port of src/cnvscope/qc.py (qc_matrices + compute_qc): same
// definitions, written one cell at a time; tests check it against the Python.

import type { Semantics } from "../api/client";
import type { Column } from "./dataset";

/** numpy's round: halves go to the even neighbour. */
function roundEven(x: number): number {
  const r = Math.round(x);
  return Math.abs(x % 1) === 0.5 && r % 2 !== 0 ? r - 1 : r;
}

/** numpy's nanmedian: NaN when no finite value. */
function nanmedian(v: number[]): number {
  const s = v.filter((x) => !Number.isNaN(x)).sort((a, b) => a - b);
  const n = s.length;
  if (!n) return NaN;
  return n % 2 ? s[(n - 1) / 2] : (s[n / 2 - 1] + s[n / 2]) / 2;
}

function gini(r: ArrayLike<number>): number {
  const s = Array.from(r).filter((x) => Number.isFinite(x) && x >= 0).sort((a, b) => a - b);
  const n = s.length;
  if (!n) return NaN;
  let total = 0, weighted = 0;
  s.forEach((x, k) => { total += x; weighted += x * (k + 1); });
  return total === 0 ? 0 : (2 * weighted) / (n * total) - (n + 1) / n;
}

function mode(p: ArrayLike<number>): number {
  const c = new Map<number, number>();
  for (let i = 0; i < p.length; i++) if (Number.isFinite(p[i])) c.set(p[i], (c.get(p[i]) ?? 0) + 1);
  let best = NaN, bc = 0;
  for (const [v, k] of c) if (k > bc || (k === bc && v < best)) { best = v; bc = k; }
  return best;
}

/** One cell's metrics; p, l, r are its rows (or null when the layer is missing). */
export function cellMetrics(
  chrs: string[], widthMb: number[],
  p: ArrayLike<number> | null, l: ArrayLike<number> | null, r: ArrayLike<number> | null,
  baseline = 2,
): Record<string, number> {
  const out: Record<string, number> = {};
  const nb = chrs.length;
  const same = (j: number) => j > 0 && chrs[j] === chrs[j - 1];
  if (r) {
    let tot = 0, n = 0;
    for (let j = 0; j < nb; j++) if (Number.isFinite(r[j])) { tot += r[j]; n++; }
    out.n_reads_total = tot;
    out.reads_per_bin = n > 0 ? tot / n : NaN;
    out.gini = gini(r);
  }
  if (!p) return out;

  const fin = (j: number) => Number.isFinite(p[j]);
  let sum = 0, nf = 0;
  for (let j = 0; j < nb; j++) if (fin(j)) { sum += p[j]; nf++; }
  const psi = nf > 0 ? sum / nf : NaN;
  out.ploidy = psi;
  let cna = 0, bp = 0;
  for (let j = 0; j < nb; j++) {
    if (fin(j) && roundEven(p[j]) !== baseline) cna++;
    if (j > 0 && fin(j) && fin(j - 1) && p[j] !== p[j - 1] && same(j)) bp++;
  }
  out.bins_with_cna = cna;
  out.breakpoints = bp;
  out.state_mode = mode(p);
  if (!l) return out;

  const okCell = Number.isFinite(psi) && psi > 0;
  const keys = ["mapd", "median_residual", "spikiness", "distance_integer_logR", "MBRSM_dispersion", "autocorr", "MSRSI_non_integerness"];
  if (!okCell) { for (const k of keys) out[k] = NaN; return out; }

  const d1: number[] = [], res: number[] = [], d2: number[] = [];
  for (let j = 0; j < nb; j++) {
    if (j > 0 && same(j)) d1.push(Math.abs(l[j] - l[j - 1]));
    res.push(Math.abs(psi * 2 ** l[j] - p[j]));
    if (j >= 2 && same(j - 1) && same(j)) d2.push(Math.abs(l[j] - 2 * l[j - 1] + l[j - 2]));
  }
  out.mapd = nanmedian(d1);
  out.median_residual = nanmedian(res);
  out.spikiness = nanmedian(d2);

  // segments: runs of one state within a chromosome (a missing state starts a run each bin)
  const seg = new Float64Array(nb).fill(NaN);
  let dsum = 0, dw = 0;
  for (let a = 0; a < nb; ) {
    let b = a + 1;
    while (b < nb && same(b) && p[b] === p[b - 1]) b++;
    let cnt = 0, sl = 0, w = 0, anyFin = false;
    for (let j = a; j < b; j++) {
      if (Number.isFinite(l[j])) { cnt++; sl += l[j]; }
      if (fin(j)) anyFin = true;
      w += widthMb[j];
    }
    if (cnt > 0 && anyFin) {
      const mean = sl / cnt;
      seg.fill(mean, a, b);
      const cn = psi * 2 ** mean;
      dsum += w * (cn - roundEven(cn)) ** 2;
      dw += w;
    }
    a = b;
  }
  out.distance_integer_logR = dw > 0 ? dsum / dw : NaN;

  const resid = Array.from({ length: nb }, (_, j) => l[j] - seg[j]);
  const med = nanmedian(resid);
  out.MBRSM_dispersion = 1.4826 * nanmedian(resid.map((x) => Math.abs(x - med)));

  // lag-1 autocorrelation of the residuals per chromosome, weighted by its bins
  let acs = 0, aw = 0;
  for (let a = 0; a < nb; ) {
    let b = a + 1;
    while (b < nb && same(b)) b++;
    if (b - a >= 3) {
      let k = 0, mx = 0, my = 0;
      for (let j = a; j < b - 1; j++) if (Number.isFinite(resid[j]) && Number.isFinite(resid[j + 1])) { k++; mx += resid[j]; my += resid[j + 1]; }
      if (k >= 2) {
        mx /= k; my /= k;
        let sxx = 0, syy = 0, sxy = 0;
        for (let j = a; j < b - 1; j++) {
          if (!(Number.isFinite(resid[j]) && Number.isFinite(resid[j + 1]))) continue;
          const dx = resid[j] - mx, dy = resid[j + 1] - my;
          sxx += dx * dx; syy += dy * dy; sxy += dx * dy;
        }
        const c = sxy / Math.sqrt(sxx * syy);
        if (sxx > 0 && syy > 0 && Number.isFinite(c)) { acs += c * (b - a); aw += b - a; }
      }
    }
    a = b;
  }
  out.autocorr = aw > 0 ? acs / aw : NaN;

  // distance of segment logR from the logR its integer state predicts, by (chromosome, state)
  const groups = new Map<string, number[]>();
  let ci = -1;
  for (let j = 0; j < nb; j++) {
    if (!same(j)) ci++;
    if (!(fin(j) && p[j] > 0 && Number.isFinite(seg[j]))) continue;
    const key = `${ci}|${p[j]}`;
    const v = seg[j] - Math.log2(p[j] / psi);
    (groups.get(key) ?? groups.set(key, []).get(key)!).push(v);
  }
  const meds = [...groups.values()].map(nanmedian);
  const mm = nanmedian(meds);
  out.MSRSI_non_integerness = meds.length ? 1.4826 * nanmedian(meds.map((m) => Math.abs(m - mm))) : NaN;
  return out;
}

const pick = (sem: Record<string, Semantics>, want: Semantics, prefer: string[] = []) =>
  prefer.find((n) => sem[n] === want) ?? Object.keys(sem).find((n) => sem[n] === want && !/_(lower|upper)$/.test(n)) ?? null;

/**
 * Add the QC metrics as cell columns (like compute_qc): values already there
 * are kept, metrics no cell has are left out.
 */
export function addQC(
  obs: Record<string, Column>, n: number, chrs: string[], widthMb: number[],
  layers: Record<string, Float32Array>, sem: Record<string, Semantics>,
): void {
  const nb = chrs.length;
  const P = pick(sem, "integer_cn", ["total_cn"]), R = pick(sem, "count", ["read_counts"]);
  const Lr = pick(sem, "continuous_ratio", ["logr"]), CS = pick(sem, "continuous_cn");
  if (!P && !R) return;
  const row = (name: string | null, i: number) => (name ? layers[name].subarray(i * nb, (i + 1) * nb) : null);
  const cols: Record<string, number[]> = {};
  for (let i = 0; i < n; i++) {
    const p = row(P, i);
    let l: ArrayLike<number> | null = row(Lr, i);
    if (!l && CS && p) {
      let s = 0, k = 0;
      for (let j = 0; j < nb; j++) if (Number.isFinite(p[j])) { s += p[j]; k++; }
      const psi = k ? s / k : NaN, cs = row(CS, i)!;
      l = Float64Array.from(cs, (v) => (v > 0 ? Math.log2(v / psi) : NaN));
    }
    const m = cellMetrics(chrs, widthMb, p, l, row(R, i));
    for (const [k, v] of Object.entries(m)) (cols[k] ??= new Array(n).fill(NaN))[i] = v;
  }
  for (const [k, col] of Object.entries(cols)) {
    if (!col.some(Number.isFinite)) continue;
    const have = obs[k];
    obs[k] = col.map((v, i) => (have && typeof have[i] === "number" && Number.isFinite(have[i]) ? (have[i] as number) : v));
  }
}
