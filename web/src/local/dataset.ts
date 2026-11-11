// A dataset held in the browser (no cnvscope server): the same Meta the server
// sends (src/cnvscope/server.py: build_meta) plus its layers in memory.

import type { ChrInfo, Dendrogram, Meta, Semantics } from "../api/client";

export type Column = (string | number | boolean | null)[];

export interface LocalData {
  meta: Meta;
  /** cells × bins, row-major, rows in dataset order */
  layers: Record<string, Float32Array>;
}

export interface Parts {
  name: string;
  tool: string;
  modality: Meta["modality"];
  genome: string | null;
  cells: string[];
  /** cell columns, each one value per cell */
  obs: Record<string, Column>;
  bins: { chr: string[]; start: number[]; end: number[]; start_cum?: number[]; end_cum?: number[] };
  layers: Record<string, Float32Array>;
  semantics: Record<string, Semantics>;
  default_layer: string;
  dendrograms?: Record<string, Dendrogram>;
}

/** Cumulative genome coordinates, chromosomes laid end to end in bin order. */
function cumulative(chr: string[], start: number[], end: number[]): [number[], number[]] {
  const sc: number[] = [], ec: number[] = [];
  let offset = 0, prevChr = "", chrEnd = 0;
  for (let i = 0; i < chr.length; i++) {
    if (chr[i] !== prevChr) { offset += chrEnd; chrEnd = 0; prevChr = chr[i]; }
    sc.push(offset + start[i]); ec.push(offset + end[i]);
    chrEnd = Math.max(chrEnd, end[i]);
  }
  return [sc, ec];
}

function chrInfo(chr: string[], sc: number[], ec: number[]): ChrInfo[] {
  const out: ChrInfo[] = [];
  for (let i = 0; i < chr.length; i++) {
    const last = out[out.length - 1];
    if (last && last.chr === chr[i]) { last.bin_end = i + 1; last.end_cum = ec[i]; }
    else out.push({ chr: chr[i], bin_start: i, bin_end: i + 1, start_cum: sc[i], end_cum: ec[i] });
  }
  return out;
}

function kind(col: Column): Meta["obs_kinds"][string] {
  const v = col.find((x) => x != null);
  return typeof v === "boolean" ? "boolean" : typeof v === "number" || v === undefined ? "numeric" : "categorical";
}

export function assemble(p: Parts): LocalData {
  const { chr, start, end } = p.bins;
  const [sc, ec] = p.bins.start_cum && p.bins.end_cum ? [p.bins.start_cum, p.bins.end_cum] : cumulative(chr, start, end);
  const names = Object.keys(p.layers);
  const intervals: Record<string, [string, string]> = {};
  for (const n of names) if (p.layers[`${n}_lower`] && p.layers[`${n}_upper`]) intervals[n] = [`${n}_lower`, `${n}_upper`];
  const bounds = new Set(Object.values(intervals).flat());
  const obs: Meta["obs"] = {}, kinds: Meta["obs_kinds"] = {};
  for (const [k, col] of Object.entries(p.obs)) {
    obs[k] = col.map((v) => (typeof v === "number" && !Number.isFinite(v) ? null : v));
    kinds[k] = kind(obs[k]);
  }
  const meta: Meta = {
    version: "browser",
    name: p.name,
    tool: p.tool,
    modality: p.modality,
    genome: p.genome,
    n_cells: p.cells.length,
    n_bins: chr.length,
    cells: p.cells,
    obs,
    obs_kinds: kinds,
    bins: { chr, start, end, start_cum: sc, end_cum: ec },
    chr_info: chrInfo(chr, sc, ec),
    layers: names.filter((n) => !bounds.has(n)).map((n) => ({ name: n, semantics: p.semantics[n] })),
    intervals,
    default_layer: p.default_layer,
    dendrograms: p.dendrograms ?? {},
  };
  return { meta, layers: p.layers };
}
