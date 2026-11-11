// Per-cell QC metrics: labels, which side is bad, and robust summaries (ported
// from the ASCAT.sc Explorer's qc.js). Metrics are the numeric cell columns the
// server sends; known ones get ASCAT.sc's getQCs() labels, others a generic one.

import { quantile } from "./canvas";

/** Side on which a value is bad; "none": a description of the cell (e.g. its CNV burden), not a quality measure. */
export type Dir = "upper" | "lower" | "both" | "none";
export interface MetricInfo {
  key: string;
  label: string;
  group: string;
  dir: Dir;
  log: boolean;
  help: string;
}
export interface Metric extends MetricInfo {
  n: number;
  median: number;
  mad: number | null;
  lo: number;
  hi: number;
  min: number;
  max: number;
  out: number;
}

// [key, label, group, direction of a bad value, log scale, help]
const META: [string, string, string, Dir, boolean, string][] = [
  ["n_reads_total", "Total reads", "Coverage", "lower", true, "Reads in the cell's bins. Few reads give noisy profiles."],
  ["reads_per_bin", "Reads per bin", "Coverage", "lower", true, "A cell's mean reads per bin."],
  ["coverage", "Coverage", "Coverage", "lower", false, "The tool's coverage of the cell (ASCAT.sc: allSolutions)."],
  ["gini", "Gini coefficient", "Coverage", "upper", false, "Unevenness of the reads across bins: 0 is uniform, towards 1 reads pile up in few bins."],
  ["mapd", "MAPD", "Noise", "upper", false, "Median absolute difference of neighbouring bins' logR, within chromosomes."],
  ["MAPD_gc_corrected", "MAPD (getQCs)", "Noise", "upper", false, "MAPD as getQCs() computes it."],
  ["spikiness", "Spikiness", "Noise", "upper", false, "Median second difference of logR: isolated outlier bins. Lower is smoother."],
  ["MBRSM_dispersion", "Segment dispersion", "Noise", "upper", false, "Scaled MAD of the bins' logR around their segment."],
  ["autocorr", "Autocorrelation", "Noise", "upper", false, "Lag-1 autocorrelation of segment residuals; high values suggest waves (GC, replication timing)."],
  ["median_residual", "Median residual", "Fit", "upper", false, "Median distance of the bins' copy number from the fitted integer profile."],
  ["distance_integer_logR", "LogR integer distance", "Fit", "upper", false, "Distance of segment copy numbers from the nearest integer. Lower is a better fit."],
  ["distance_integer_BAF", "BAF integer distance", "Fit", "upper", false, "Distance of segment BAFs from the nearest integer allelic state. Lower is cleaner."],
  ["MSRSI_non_integerness", "Non-integerness", "Fit", "upper", false, "Scaled MAD of segment logR from the logR its integer state predicts."],
  ["breakpoints", "Breakpoints", "Fit", "upper", false, "Copy-number state changes; very many usually means noise was segmented."],
  ["bins_with_cna", "Bins with CNA", "Solution", "upper", false, "Bins whose copy number is not the baseline (after the minimum segment size)."],
  ["state_mode", "Modal CN state", "Solution", "both", false, "Most frequent copy-number state; usually 2."],
  ["ploidy", "Ploidy", "Solution", "both", false, "Mean copy number over the bins."],
  ["purity", "Purity", "Solution", "both", false, "Fitted purity; single cells are expected near 1."],
  ["cnv_burden", "CNV burden", "Karyotype", "none", false, "Fraction of the genome's windows called as gain or loss (epiAneufinder's cnv_burden_per_cell)."],
  ["gain_fraction", "Gained fraction", "Karyotype", "none", false, "Fraction of the windows called as gain."],
  ["loss_fraction", "Lost fraction", "Karyotype", "none", false, "Fraction of the windows called as loss."],
];
export const KNOWN: Record<string, MetricInfo> = Object.fromEntries(
  META.map(([key, label, group, dir, log, help]) => [key, { key, label, group, dir, log, help }]),
);
export const INT_METRICS = new Set(["breakpoints", "state_mode", "bins_with_cna", "n_reads_total"]);

export function infoFor(key: string): MetricInfo {
  return KNOWN[key] ?? { key, label: key.replace(/_/g, " "), group: "Other", dir: "both", log: false, help: "A per-cell value from the dataset." };
}

/** Robust summary of a metric over all cells, with median ± 3 MAD outliers on its bad side. */
export function summarise(m: MetricInfo, values: ArrayLike<number>): Metric {
  const vs = Array.from(values);
  const sc = vs.map((v) => (m.log ? (v > 0 ? Math.log10(v) : NaN) : v));
  const ok = sc.filter((v) => !Number.isNaN(v)).sort((a, b) => a - b);
  const med = quantile(ok, 0.5), mad = 1.4826 * quantile(ok.map((v) => Math.abs(v - med)).sort((a, b) => a - b), 0.5);
  const lo = med - 3 * mad, hi = med + 3 * mad, back = (v: number) => (m.log ? Math.pow(10, v) : v);
  let out = 0;
  sc.forEach((v) => {
    if (Number.isNaN(v) || !(mad > 0) || m.dir === "none") return;
    if ((m.dir !== "lower" && v > hi) || (m.dir !== "upper" && v < lo)) out++;
  });
  const fin = vs.filter((v) => !Number.isNaN(v));
  let min = Infinity, max = -Infinity;
  for (const v of fin) { if (v < min) min = v; if (v > max) max = v; }
  return { ...m, n: fin.length, median: back(med), mad: m.log ? null : mad, lo: back(lo), hi: back(hi), min, max, out };
}

/** Summaries of the metrics in `vals`: known ones first, in their usual order. */
export function metricSummaries(vals: Record<string, Float64Array>): Metric[] {
  const keys = Object.keys(vals).filter((k) => vals[k].some((v) => !Number.isNaN(v)));
  const order = (k: string) => { const i = META.findIndex(([x]) => x === k); return i < 0 ? 1e3 : i; };
  keys.sort((a, b) => order(a) - order(b) || a.localeCompare(b));
  return keys.map((k) => summarise(infoFor(k), vals[k]));
}

/** Whether v lies beyond median ± 3 MAD on the metric's bad side. */
export function isOutlier(m: Metric, v: number): boolean {
  if (v == null || Number.isNaN(v) || m.dir === "none") return false;
  if (m.log && !(v > 0)) return m.dir !== "upper";
  return (m.dir !== "lower" && v > m.hi) || (m.dir !== "upper" && v < m.lo);
}

export const fmtMetric = (key: string, v: number | null | undefined): string => {
  if (v == null || !Number.isFinite(v)) return "–";
  if (key === "n_reads_total") return v < 1e3 ? String(Math.round(v)) : v < 1e6 ? (v / 1e3).toFixed(v < 1e4 ? 2 : 1) + "k" : (v / 1e6).toFixed(v < 1e7 ? 2 : 1) + "M";
  if (INT_METRICS.has(key)) return String(Math.round(v * 10) / 10);
  const a = Math.abs(v);
  return a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) : v.toFixed(3);
};

/* ---------- sidebar QC thresholds: a cell beyond any active one is flagged ---------- */
export interface Threshold {
  key: string;
  metric: string;
  label: string;
  dir: "max" | "min";
  reason: string;
  tile: string;
  hint: string;
  f: (v: number) => string;
  int?: boolean;
  /** default for a new dataset */
  initial: (sorted: number[]) => number | null;
}

const robust = (v: number[]) => {
  if (!v.length) return null;
  const m = quantile(v, 0.5), mad = 1.4826 * quantile(v.map((x) => Math.abs(x - m)).sort((a, b) => a - b), 0.5);
  return mad > 0 ? Math.round((m + 2 * mad) * 1000) / 1000 : null;
};
const pct = (p: number, f: (v: number) => number) => (v: number[]) => (v.length > 1 && v[0] < v[v.length - 1] ? f(quantile(v, p)) : null);

// Defaults: median + 2 MAD for the residual and MAPD, the 10th percentile of
// coverage, the 90th of CNA bins and Gini.
export const THRESH: Threshold[] = [
  { key: "residual", metric: "median_residual", label: "Max median residual", dir: "max", reason: "high residual", tile: "High residual", hint: "distance of the bins from the integer fit", f: (v) => v.toFixed(2), initial: robust },
  { key: "mapd", metric: "mapd", label: "Max MAPD", dir: "max", reason: "noisy", tile: "Noisy", hint: "noise of the cell's logR", f: (v) => v.toFixed(3), initial: robust },
  { key: "coverage", metric: "coverage", label: "Min coverage", dir: "min", reason: "low coverage", tile: "Low coverage", hint: "the tool's coverage of the cell", f: (v) => Math.round(v) + "×", int: true, initial: pct(0.1, Math.floor) },
  { key: "bins_with_cna", metric: "bins_with_cna", label: "Max bins with CNA", dir: "max", reason: "many CNA bins", tile: "Many CNA bins", hint: "bins whose copy number is not the baseline", f: (v) => String(Math.round(v)), int: true, initial: pct(0.9, Math.ceil) },
  { key: "gini", metric: "gini", label: "Max Gini", dir: "max", reason: "uneven coverage", tile: "Uneven coverage", hint: "unevenness of the reads across bins", f: (v) => v.toFixed(3), initial: pct(0.9, (v) => Math.round(v * 1000) / 1000) },
];
