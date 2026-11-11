// Read a cnvscope .h5ad (AnnData, HDF5) in the browser with h5wasm.
// Mirrors CNVDataset.read + build_meta: obs/var dataframes, layers, uns.

import type { Dendrogram, Meta, Semantics } from "../api/client";
import { assemble, type Column, type LocalData } from "./dataset";

type H5 = typeof import("h5wasm");
// eslint-disable-next-line @typescript-eslint/no-explicit-any
type Node = any;

let lib: Promise<H5> | null = null;
/** h5wasm (about 1 MB of WebAssembly), loaded the first time a file is opened. */
function h5(): Promise<H5> {
  lib ??= import("h5wasm").then(async (m) => { await m.ready; return m; });
  return lib;
}

const isGroup = (n: Node) => typeof n?.keys === "function";
const enc = (n: Node): string => String(n?.attrs?.["encoding-type"]?.value ?? "");
const asList = (v: unknown): unknown[] => (v == null ? [] : Array.isArray(v) || ArrayBuffer.isView(v) ? Array.from(v as ArrayLike<unknown>) : [v]);
const num = (v: unknown) => (typeof v === "bigint" ? Number(v) : (v as number));

/** One element of an AnnData file as plain JS values. */
function read(n: Node): unknown {
  if (!n) return null;
  if (!isGroup(n)) {
    const v = n.value;
    if (v instanceof BigInt64Array || v instanceof BigUint64Array) return Array.from(v, Number);
    return v;
  }
  const e = enc(n);
  if (e === "categorical") {
    const cats = asList(read(n.get("categories"))).map(String);
    return asList(read(n.get("codes"))).map((c) => (num(c) < 0 ? null : cats[num(c)]));
  }
  if (n.keys().includes("values") && n.keys().includes("mask")) {
    const vals = asList(read(n.get("values"))), mask = asList(read(n.get("mask")));
    return vals.map((v, i) => (mask[i] ? null : typeof v === "bigint" ? Number(v) : v));
  }
  const out: Record<string, unknown> = {};
  for (const k of n.keys()) out[k] = read(n.get(k));
  return out;
}

function column(v: unknown): Column {
  return asList(v).map((x) => {
    if (x == null) return null;
    if (typeof x === "bigint") return Number(x);
    if (typeof x === "number" || typeof x === "boolean") return x;
    return String(x);
  });
}

function frame(g: Node): { index: string[]; cols: Record<string, Column> } {
  const idxName = String(g.attrs["_index"]?.value ?? "_index");
  const order = asList(g.attrs["column-order"]?.value).map(String);
  const cols: Record<string, Column> = {};
  for (const k of order) cols[k] = column(read(g.get(k)));
  return { index: asList(read(g.get(idxName))).map(String), cols };
}

export async function readH5ad(bytes: Uint8Array, name: string): Promise<LocalData> {
  const m = await h5();
  const path = `/cnvscope-${Date.now()}.h5ad`;
  m.FS!.writeFile(path, bytes);
  const f = new m.File(path, "r");
  try {
    const obs = frame(f.get("obs")), varF = frame(f.get("var"));
    const uns = (read(f.get("uns")) ?? {}) as Record<string, unknown>;
    const semantics = (uns.semantics ?? {}) as Record<string, Semantics>;
    const lg: Node = f.get("layers");
    const layers: Record<string, Float32Array> = {};
    for (const k of lg?.keys() ?? []) {
      const d = lg.get(k);
      if (isGroup(d)) throw new Error(`layer ${k} is stored sparse; cnvscope expects dense layers`);
      const v = d.value;
      layers[k] = v instanceof Float32Array ? v : Float32Array.from(v as ArrayLike<number>);
    }
    const dendrograms: Record<string, Dendrogram> = {};
    for (const [k, d] of Object.entries((uns.dendrograms ?? {}) as Record<string, Record<string, unknown>>)) {
      const merge = asList(d.merge).map(num);
      dendrograms[k] = {
        merge: Array.from({ length: merge.length / 2 }, (_, i) => [merge[2 * i], merge[2 * i + 1]] as [number, number]),
        height: asList(d.height).map(num),
        order: asList(d.order).map(num),
        layers: asList(d.layers).map(String),
        method: String(d.method ?? ""),
      };
    }
    const str = (k: string) => (uns[k] == null ? null : String(uns[k]));
    const ints = (k: string) => (varF.cols[k] ?? []).map((v) => Number(v));
    return assemble({
      name: name || str("source") || str("tool") || "dataset",
      tool: str("tool") ?? "",
      modality: (str("modality") ?? "scDNA") as Meta["modality"],
      genome: str("genome"),
      cells: obs.index,
      obs: obs.cols,
      bins: {
        chr: (varF.cols.chr ?? []).map(String),
        start: ints("start"),
        end: ints("end"),
        start_cum: varF.cols.start_cum ? ints("start_cum") : undefined,
        end_cum: varF.cols.end_cum ? ints("end_cum") : undefined,
      },
      layers,
      semantics,
      default_layer: str("default_layer") ?? Object.keys(layers)[0],
      dendrograms,
    });
  } finally {
    f.close();
    m.FS!.unlink(path);
  }
}
