// QC filters: any per-cell metric can have a minimum, a maximum or both. A cell outside any
// active filter is flagged. Defaults come from THRESH (median + 2 MAD for the residual and
// MAPD, percentiles for coverage, CNA bins and Gini), the ASCAT.sc Explorer's thresholds.

import type { Meta } from "../api/client";
import { THRESH, fmtMetric, infoFor } from "./metrics";

export interface Filter {
  min: number | null;
  max: number | null;
}
export type Filters = Record<string, Filter>;

export const isActive = (f: Filter | undefined): f is Filter => !!f && (f.min != null || f.max != null);

/** Finite min and max of each metric, for the ones that vary. */
export function metricRanges(vals: Record<string, ArrayLike<number>>): Record<string, [number, number]> {
  const out: Record<string, [number, number]> = {};
  for (const [k, v] of Object.entries(vals)) {
    let lo = Infinity, hi = -Infinity;
    for (let i = 0; i < v.length; i++) { const x = v[i]; if (Number.isFinite(x)) { if (x < lo) lo = x; if (x > hi) hi = x; } }
    if (lo < hi) out[k] = [lo, hi];
  }
  return out;
}

/** The ASCAT.sc Explorer's default thresholds, for the metrics this dataset has. */
export function defaultFilters(vals: Record<string, ArrayLike<number>>): Filters {
  const ranges = metricRanges(vals);
  const out: Filters = {};
  for (const t of THRESH) {
    if (!ranges[t.metric]) continue;
    const v = Array.from(vals[t.metric]).filter(Number.isFinite).sort((a, b) => a - b);
    const x = t.initial(v);
    // skip a default that cannot flag any cell (at or beyond the data's own end)
    if (x == null || (t.dir === "max" ? x >= ranges[t.metric][1] : x <= ranges[t.metric][0])) continue;
    out[t.metric] = t.dir === "max" ? { min: null, max: x } : { min: x, max: null };
  }
  return out;
}

/** Why a cell fails one side of a filter, e.g. "noisy" or "low ploidy". */
export function reasonFor(metric: string, side: "min" | "max"): string {
  const t = THRESH.find((t) => t.metric === metric && (t.dir === "max") === (side === "max"));
  return t ? t.reason : `${side === "max" ? "high" : "low"} ${infoFor(metric).label.toLowerCase()}`;
}

export interface Evaluation {
  passing: boolean[];
  reasons: string[];
  /** cells flagged by each metric's filter */
  failCounts: Record<string, number>;
  nPass: number;
}

export function evaluate(n: number, vals: Record<string, ArrayLike<number>>, filters: Filters): Evaluation {
  const passing = new Array<boolean>(n).fill(true), reasons = new Array<string>(n).fill("");
  const failCounts: Record<string, number> = {};
  for (const [metric, f] of Object.entries(filters)) {
    const v = vals[metric];
    if (!isActive(f) || !v) continue;
    failCounts[metric] = 0;
    for (let i = 0; i < n; i++) {
      const x = v[i];
      if (!Number.isFinite(x)) continue;
      const side = f.min != null && x < f.min ? "min" : f.max != null && x > f.max ? "max" : null;
      if (!side) continue;
      passing[i] = false;
      const r = reasonFor(metric, side);
      reasons[i] = reasons[i] ? reasons[i] + ", " + r : r;
      failCounts[metric]++;
    }
  }
  return { passing, reasons, failCounts, nPass: passing.filter(Boolean).length };
}

/** Whether a value fails a metric's filter. */
export function fails(f: Filter | undefined, v: number): boolean {
  return isActive(f) && Number.isFinite(v) && ((f.min != null && v < f.min) || (f.max != null && v > f.max));
}

/** "≤ 0.29", "≥ 3", "0.2 – 0.5" */
export function describe(metric: string, f: Filter): string {
  const m = (v: number) => fmtMetric(metric, v);
  if (f.min != null && f.max != null) return `${m(f.min)} – ${m(f.max)}`;
  if (f.max != null) return `≤ ${m(f.max)}`;
  return `≥ ${m(f.min!)}`;
}

/** A tab-separated table of every cell: pass/fail, failed filters, all metrics. */
export function cellsTable(meta: Meta, vals: Record<string, ArrayLike<number>>, filters: Filters, ev: Evaluation): string {
  const metrics = Object.keys(vals).sort();
  const annots = Object.keys(meta.obs_kinds).filter((k) => meta.obs_kinds[k] !== "numeric");
  const active = Object.entries(filters).filter(([, f]) => isActive(f));
  const lines = [
    `# cnvscope ${meta.version} · ${meta.name} · ${meta.tool} · ${new Date().toISOString().slice(0, 19)}Z`,
    `# ${ev.nPass} of ${meta.n_cells} cells pass the filters`,
    ...(active.length ? active.map(([k, f]) => `# filter: ${k} ${describe(k, f)}`) : ["# filter: none"]),
    ["cell", "passes_qc", "failed", ...annots, ...metrics].join("\t"),
  ];
  const cell = (x: unknown) => (x == null || (typeof x === "number" && !Number.isFinite(x)) ? "" : String(x).replace(/[\t\n]/g, " "));
  for (let i = 0; i < meta.n_cells; i++) {
    lines.push([meta.cells[i], ev.passing[i] ? "TRUE" : "FALSE", ev.reasons[i], ...annots.map((a) => cell(meta.obs[a][i])), ...metrics.map((m) => cell(vals[m][i]))].map(cell).join("\t"));
  }
  return lines.join("\n") + "\n";
}

export function downloadText(text: string, filename: string, type = "text/tab-separated-values") {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}
