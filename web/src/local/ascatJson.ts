// ASCAT.sc web export (ascat_data.json[.gz], from rds_to_web() or export_cnvscope())
// read in the browser. A port of src/cnvscope/readers/ascatsc_json.py: from_dict.

import type { Dendrogram, Semantics } from "../api/client";
import { layoutOrder } from "./hclust";
import { assemble, type Column, type LocalData } from "./dataset";
import { addQC } from "./qc";

type Vals = (number | null)[];
// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Json = any;

/** The file's text, gunzipped when it is gzip. */
export async function fileText(bytes: Uint8Array): Promise<string> {
  if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
    const s = new Blob([bytes as BlobPart]).stream().pipeThrough(new DecompressionStream("gzip"));
    return new Response(s).text();
  }
  return new TextDecoder().decode(bytes);
}

function matrix(d: Record<string, Vals> | null | undefined, cells: string[], nb: number): Float32Array | null {
  if (!d || !Object.keys(d).length) return null;
  const m = new Float32Array(cells.length * nb).fill(NaN);
  cells.forEach((c, i) => {
    const v = d[c];
    if (v) for (let j = 0; j < nb; j++) m[i * nb + j] = v[j] == null ? NaN : (v[j] as number);
  });
  return m;
}

export function fromDict(d: Json, name: string): LocalData {
  if (!d?.bins || !d?.profiles) throw new Error("not an ASCAT.sc web export: missing 'bins' or 'profiles'");
  const cells = Object.keys(d.profiles);
  const b = d.bins, nb = b.chr.length;
  const layers: Record<string, Float32Array> = { total_cn: matrix(d.profiles, cells, nb)! };
  const sem: Record<string, Semantics> = { total_cn: "integer_cn" };
  if (d.nMajor && d.nMinor && Object.keys(d.nMajor).length) {
    layers.major_cn = matrix(d.nMajor, cells, nb)!;
    layers.minor_cn = matrix(d.nMinor, cells, nb)!;
    sem.major_cn = sem.minor_cn = "allele_cn";
    const ci = d.ci ?? {};
    if (Object.keys(ci).length)
      for (const [allele, base] of [["nMajor", "major_cn"], ["nMinor", "minor_cn"]])
        for (const side of ["lower", "upper"]) {
          const key = `${allele}_${side}`, sub: Record<string, Vals> = {};
          for (const [c, v] of Object.entries(ci as Record<string, Record<string, Vals>>)) if (v?.[key] != null) sub[c] = v[key];
          const m = matrix(sub, cells, nb);
          if (m) { layers[`${base}_${side}`] = m; sem[`${base}_${side}`] = "allele_cn"; }
        }
  }
  const raw = matrix(d.raw, cells, nb);
  if (raw) {
    let min = Infinity;
    for (let i = 0; i < Math.min(20, cells.length) * nb; i++) if (raw[i] < min) min = raw[i];
    if (min < 0) { layers.logr = raw; sem.logr = "continuous_ratio"; } else { layers.cn_scale = raw; sem.cn_scale = "continuous_cn"; }
  }
  const reads = matrix(d.reads, cells, nb);
  if (reads) { layers.read_counts = reads; sem.read_counts = "count"; }

  const obs: Record<string, Column> = {};
  for (const src of ["quality", "qc_metrics"])
    for (const [c, q] of Object.entries((d[src] ?? {}) as Record<string, Record<string, unknown>>)) {
      const i = cells.indexOf(c);
      if (i < 0 || !q || typeof q !== "object") continue;
      for (const [k, v] of Object.entries(q))
        if (typeof v === "number") (obs[k] ??= new Array(cells.length).fill(null))[i] = v;
    }
  if (d.cell_types && Object.keys(d.cell_types).length) obs.cell_type = cells.map((c) => d.cell_types[c] ?? null);

  const dendrograms: Record<string, Dendrogram> = {};
  const dend = d.dendrogram;
  if (dend?.merge?.length === cells.length - 1) {
    const merge = dend.merge as [number, number][];
    let order = layoutOrder(merge);
    const co = d.clustering_order;
    if (co?.length === cells.length) {
      if (co.every((x: unknown) => typeof x === "string")) {
        const idx = new Map(cells.map((c, i) => [c, i]));
        if (co.every((c: string) => idx.has(c))) order = co.map((c: string) => idx.get(c)!);
      } else if (co.every((x: unknown) => typeof x === "number")) {
        const o = (co as number[]).map((x) => x - (Math.min(...co) === 1 ? 1 : 0));
        if ([...o].sort((a, b) => a - b).every((v, i) => v === i)) order = o;
      }
    }
    dendrograms.total = { merge, height: dend.height, order, layers: ["total_cn"], method: "ward.D2, Manhattan (from the export)" };
  }

  const width = (b.end as number[]).map((e, j) => (e - b.start[j]) / 1e6);
  addQC(obs, cells.length, (b.chr as unknown[]).map(String), width, layers, sem);
  return assemble({
    name, tool: "ASCAT.sc", modality: "scDNA", genome: "hg38", cells, obs,
    bins: { chr: b.chr.map(String), start: b.start, end: b.end, start_cum: b.start_cum, end_cum: b.end_cum },
    layers, semantics: sem, default_layer: "total_cn", dendrograms,
  });
}

export async function readAscatJson(bytes: Uint8Array, name: string): Promise<LocalData> {
  return fromDict(JSON.parse(await fileText(bytes)), name);
}
