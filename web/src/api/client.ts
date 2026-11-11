// Client for the cnvscope API (src/cnvscope/server.py), or, with no server
// (the GitHub Pages site, or a build made with VITE_STATIC=1), the same calls
// answered in the browser from a file read locally (src/local/).

import type { LocalData } from "../local/dataset";

export type Semantics =
  | "integer_cn"
  | "allele_cn"
  | "ternary_state"
  | "continuous_ratio"
  | "continuous_cn"
  | "ratio"
  | "baf"
  | "probability"
  | "count";

export interface ChrInfo {
  chr: string;
  bin_start: number;
  bin_end: number;
  start_cum: number;
  end_cum: number;
}

export interface Dendrogram {
  /** R hclust layout over dataset row indices */
  merge: [number, number][];
  height: number[];
  /** leaf order, dataset row indices */
  order: number[];
  layers: string[];
  method: string;
}

export interface Meta {
  version: string;
  name: string;
  tool: string;
  modality: "scDNA" | "scATAC" | "scRNA" | "multiome";
  genome: string | null;
  n_cells: number;
  n_bins: number;
  cells: string[];
  obs: Record<string, (string | number | boolean | null)[]>;
  obs_kinds: Record<string, "categorical" | "numeric" | "boolean">;
  bins: { chr: string[]; start: number[]; end: number[]; start_cum: number[]; end_cum: number[] };
  chr_info: ChrInfo[];
  layers: { name: string; semantics: Semantics }[];
  /** layer → [lower, upper] interval layers */
  intervals: Record<string, [string, string]>;
  default_layer: string;
  dendrograms: Record<string, Dendrogram>;
}

export interface CellDetail {
  cell: string;
  layers: Record<string, (number | null)[]>;
}

async function json<T>(r: Response): Promise<T> {
  if (!r.ok) {
    let msg = `HTTP ${r.status}`;
    try {
      const body = await r.json();
      if (body?.detail) msg = typeof body.detail === "string" ? body.detail : JSON.stringify(body.detail);
    } catch {
      /* not JSON */
    }
    throw new Error(msg);
  }
  return r.json();
}

/** No cnvscope server: files are read in the browser and never leave it. */
let browserOnly = import.meta.env.VITE_STATIC === "1";
let local: LocalData | null = null;
export const isBrowserOnly = () => browserOnly;

export async function fetchMeta(): Promise<Meta | null> {
  if (!browserOnly) {
    try {
      const r = await fetch("api/meta");
      if (r.ok && (r.headers.get("content-type") ?? "").includes("json")) return (await json<{ dataset: Meta | null }>(r)).dataset;
    } catch {
      /* no server */
    }
    browserOnly = true;
  }
  return local?.meta ?? null;
}

async function readLocal(bytes: Uint8Array, fileName: string, name: string): Promise<Meta> {
  const n = fileName.toLowerCase();
  if (n.endsWith(".h5ad")) {
    const { readH5ad } = await import("../local/h5ad");
    local = await readH5ad(bytes, name);
  } else if (n.endsWith(".json") || n.endsWith(".json.gz") || n.endsWith(".gz")) {
    const { readAscatJson } = await import("../local/ascatJson");
    local = await readAscatJson(bytes, name);
  } else {
    throw new Error(
      `${fileName}: the browser opens .h5ad files and ASCAT.sc .json exports. Read other tools' outputs ` +
        "(and .rds/.Rda) with the cnvscope command, e.g. cnvscope --mode dna --tool chisel --input out/ --outdir results/, " +
        "and open the .h5ad it writes.",
    );
  }
  return local.meta;
}

export async function openFile(file: File): Promise<Meta> {
  if (browserOnly) return readLocal(new Uint8Array(await file.arrayBuffer()), file.name, file.name.split(".")[0]);
  const r = await fetch(`api/open?name=${encodeURIComponent(file.name)}`, { method: "POST", body: file });
  return (await json<{ dataset: Meta }>(r)).dataset;
}

export async function openDemo(): Promise<Meta> {
  if (browserOnly) {
    const r = await fetch("demo.h5ad");
    if (!r.ok) throw new Error(`demo dataset: HTTP ${r.status}`);
    return readLocal(new Uint8Array(await r.arrayBuffer()), "demo.h5ad", "Demo dataset");
  }
  return (await json<{ dataset: Meta }>(await fetch("api/demo", { method: "POST" }))).dataset;
}

export async function fetchCell(index: number): Promise<CellDetail> {
  if (browserOnly) {
    if (!local) throw new Error("no dataset is open");
    const nb = local.meta.n_bins, layers: CellDetail["layers"] = {};
    for (const [k, m] of Object.entries(local.layers))
      layers[k] = Array.from(m.subarray(index * nb, (index + 1) * nb), (v) => (Number.isFinite(v) ? v : null));
    return { cell: local.meta.cells[index], layers };
  }
  return json(await fetch(`api/cell/${index}`));
}

/** Bytes per request; rows are fetched in chunks of about this size. */
const CHUNK_BYTES = 16 * 1024 * 1024;

/** A whole layer, rows in dataset order, row-major. Fetched in row chunks. */
export async function fetchLayer(name: string, nCells: number, nBins: number, onProgress?: (rowsDone: number) => void): Promise<Float32Array> {
  if (browserOnly) {
    const m = local?.layers[name];
    if (!m) throw new Error(`no layer ${name}`);
    onProgress?.(nCells);
    return m;
  }
  const out = new Float32Array(nCells * nBins);
  const step = Math.max(1, Math.floor(CHUNK_BYTES / (4 * Math.max(1, nBins))));
  for (let start = 0; start < nCells; start += step) {
    const stop = Math.min(nCells, start + step);
    const r = await fetch(`api/layer/${encodeURIComponent(name)}?start=${start}&stop=${stop}`);
    if (!r.ok) throw new Error(`layer ${name}: HTTP ${r.status}`);
    out.set(new Float32Array(await r.arrayBuffer()), start * nBins);
    onProgress?.(stop);
  }
  return out;
}
