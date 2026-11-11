// Everything computed from the dataset and the view settings: QC flags, the
// cells shown and their order, the heatmap's matrix and colour scale. Computed
// once in <DerivedProvider> and read anywhere with useDerived().

import { createContext, useContext, useDeferredValue, useMemo, type ReactNode } from "react";
import type { Dendrogram, Meta, Semantics } from "../api/client";
import { asStateMatrix } from "../lib/alleles";
import { evaluate, metricRanges, type Filters } from "../lib/filters";
import { metricSummaries, type Metric } from "../lib/metrics";
import { PALETTES } from "../lib/palettes";
import { alleleScale, isCall, scaleFor, type Scale } from "../lib/scales";
import { binsAltered, filterSmallSegments, followFilter } from "../lib/segments";
import { ALLELE, useStore } from "./store";

export interface Derived {
  meta: Meta;
  sem: Record<string, Semantics>;
  /** the integer copy-number layer, if any */
  totalLayer: string | null;
  countLayer: string | null;
  /** per-bin values to plot on the copy-number scale: [layer, how] */
  binsLayer: [string, "cn" | "logr"] | null;
  hasAS: boolean;
  hasCi: boolean;
  binSizeMb: number;
  minBins: number;
  cellTypes: string[];
  /** numeric per-cell values by metric, aligned with dataset rows */
  vals: Record<string, Float64Array>;
  metAll: Metric[];
  /** min and max of every metric that varies */
  ranges: Record<string, [number, number]>;
  filters: Filters;
  passing: boolean[];
  reasons: string[];
  /** cells flagged by each metric's filter */
  failCounts: Record<string, number>;
  nPass: number;
  dendro: Dendrogram | null;
  /** dataset rows shown, in display order */
  shown: number[];
  /** dataset row → display row, or -1 when hidden */
  rowOf: Int32Array;
  heat: { matrix: Float32Array | null; scale: Scale; progress: number | null; label: string };
}

const Ctx = createContext<Derived | null>(null);
export const useDerived = () => useContext(Ctx)!;

function pickLayer(meta: Meta, sem: Semantics, prefer: string[]): string | null {
  for (const p of prefer) if (meta.layers.some((l) => l.name === p && l.semantics === sem)) return p;
  return meta.layers.find((l) => l.semantics === sem)?.name ?? null;
}

function p99(m: Float32Array): number {
  const step = Math.max(1, Math.floor(m.length / 200000));
  const v: number[] = [];
  for (let i = 0; i < m.length; i += step) if (Number.isFinite(m[i])) v.push(m[i]);
  v.sort((a, b) => a - b);
  return v.length ? v[Math.floor(0.99 * (v.length - 1))] : 1;
}

export function DerivedProvider({ meta, children }: { meta: Meta; children: ReactNode }) {
  const view = useStore((s) => s.view);
  const matrices = useStore((s) => s.matrices);
  const loading = useStore((s) => s.loading);
  const minSegmentMb = useStore((s) => s.minSegmentMb);
  const filters = useStore((s) => s.filters);
  const passOnly = useStore((s) => s.passOnly);
  const cellType = useStore((s) => s.cellType);
  const palette = PALETTES[useStore((s) => s.palette)];
  const dark = useStore((s) => s.theme) === "dark";

  const base = useMemo(() => {
    const sem = Object.fromEntries(meta.layers.map((l) => [l.name, l.semantics])) as Record<string, Semantics>;
    const totalLayer = pickLayer(meta, "integer_cn", [meta.default_layer, "total_cn"]);
    const countLayer = pickLayer(meta, "count", ["read_counts"]);
    const cnScale = pickLayer(meta, "continuous_cn", []);
    const ratio = pickLayer(meta, "continuous_ratio", ["logr"]);
    const binsLayer: [string, "cn" | "logr"] | null = cnScale ? [cnScale, "cn"] : ratio ? [ratio, "logr"] : null;
    const hasAS = sem.major_cn === "allele_cn" && sem.minor_cn === "allele_cn";
    const hasCi = hasAS && !!meta.intervals.major_cn;
    const b = meta.bins;
    const binSizeMb = b.start.length ? (b.end[0] - b.start[0]) / 1e6 || 1 : 1;
    const ct = meta.obs.cell_type;
    const cellTypes = ct && meta.obs_kinds.cell_type === "categorical"
      ? [...new Set(ct.filter((v): v is string => v != null).map(String))].sort()
      : [];
    return { sem, totalLayer, countLayer, binsLayer, hasAS, hasCi, binSizeMb, cellTypes };
  }, [meta]);

  const minBinsNow = minSegmentMb > 0 ? Math.max(1, Math.ceil(minSegmentMb / base.binSizeMb)) : 0;
  // the slider moves at once; the heavy per-cell work follows
  const minBins = useDeferredValue(minBinsNow);

  const total = base.totalLayer ? matrices[base.totalLayer] : undefined;
  const filteredTotal = useMemo(
    () => (total && minBins > 1 ? filterSmallSegments(total, meta.n_cells, meta.n_bins, minBins, 2) : total ?? null),
    [total, minBins, meta],
  );

  const vals = useMemo(() => {
    const out: Record<string, Float64Array> = {};
    for (const [k, kind] of Object.entries(meta.obs_kinds)) {
      if (kind !== "numeric") continue;
      const col = meta.obs[k];
      out[k] = Float64Array.from(col, (v) => (typeof v === "number" ? v : NaN));
    }
    if (minBins > 1 && filteredTotal) out.bins_with_cna = binsAltered(filteredTotal, meta.n_cells, meta.n_bins, 2);
    return out;
  }, [meta, minBins, filteredTotal]);

  const metAll = useMemo(() => metricSummaries(vals), [vals]);

  const ranges = useMemo(() => metricRanges(vals), [vals]);
  const activeFilters: Filters = useMemo(() => filters ?? {}, [filters]);
  const flags = useMemo(() => evaluate(meta.n_cells, vals, activeFilters), [meta, vals, activeFilters]);

  const dendro = useMemo(
    () => (view === ALLELE ? meta.dendrograms.allele : undefined) ?? meta.dendrograms.total ?? Object.values(meta.dendrograms)[0] ?? null,
    [meta, view],
  );

  const { shown, rowOf } = useMemo(() => {
    const order = dendro?.order ?? meta.cells.map((_, i) => i);
    const ct = meta.obs.cell_type;
    const shown = order.filter((i) => (!passOnly || flags.passing[i]) && (cellType === "All" || String(ct?.[i]) === cellType));
    const rowOf = new Int32Array(meta.n_cells).fill(-1);
    shown.forEach((r, d) => (rowOf[r] = d));
    return { shown, rowOf };
  }, [meta, dendro, passOnly, flags.passing, cellType]);

  const heat = useMemo(() => {
    const progress = (names: string[]) => {
      const p = names.filter((n) => n in loading).map((n) => loading[n] / meta.n_cells);
      return p.length ? Math.min(...p) : null;
    };
    if (view === ALLELE) {
      const maj = matrices.major_cn, min = matrices.minor_cn;
      if (!maj || !min) return { matrix: null, scale: alleleScale(palette, dark), progress: progress(["major_cn", "minor_cn"]), label: "major|minor" };
      let fMaj = maj, fMin = min;
      if (minBins > 1 && total && filteredTotal) {
        fMaj = followFilter(maj, total, filteredTotal, 2, 1);
        fMin = followFilter(min, total, filteredTotal, 2, 1);
      }
      return { matrix: asStateMatrix(fMaj, fMin, palette), scale: alleleScale(palette, dark), progress: null, label: "major|minor" };
    }
    const s = base.sem[view] ?? "integer_cn";
    const m = matrices[view];
    if (!m) return { matrix: null, scale: scaleFor(s, { palette, dark }), progress: progress([view]), label: view };
    const scale = scaleFor(s, { palette, dark, countMax: s === "count" ? p99(m) : undefined });
    let matrix = m;
    if (minBins > 1 && isCall(s)) matrix = view === base.totalLayer && filteredTotal ? filteredTotal : filterSmallSegments(m, meta.n_cells, meta.n_bins, minBins, scale.baseline);
    return { matrix, scale, progress: null, label: view };
  }, [view, matrices, loading, minBins, total, filteredTotal, base, meta, palette, dark]);

  const value: Derived = {
    meta,
    ...base,
    minBins,
    vals,
    metAll,
    ranges,
    filters: activeFilters,
    ...flags,
    dendro,
    shown,
    rowOf,
    heat,
  };
  return <Ctx.Provider value={value}>{children}</Ctx.Provider>;
}
