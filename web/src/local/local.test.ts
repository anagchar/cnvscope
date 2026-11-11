// The browser readers give what the Python side gives (scripts/make_web_fixtures.py).
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { describe, expect, it } from "vitest";
import type { Meta } from "../api/client";
import { readAscatJson } from "./ascatJson";
import type { LocalData } from "./dataset";
import { readH5ad } from "./h5ad";

const FIX = new URL("../../../tests/fixtures/", import.meta.url);
const bytes = (p: string) => new Uint8Array(readFileSync(new URL(p, FIX)));
const expected = (p: string): Meta & { _layer_sums: Record<string, number> } => JSON.parse(gunzipSync(readFileSync(new URL(p, FIX))).toString());

function close(a: unknown, b: unknown, k: string) {
  if (typeof b === "number" && typeof a === "number") expect(a, k).toBeCloseTo(b, Math.abs(b) > 1e3 ? 3 : 8);
  else expect(a, k).toEqual(b);
}

function same(got: LocalData, want: ReturnType<typeof expected>) {
  const m = got.meta;
  for (const k of ["name", "tool", "modality", "genome", "n_cells", "n_bins", "cells", "bins", "chr_info", "intervals", "default_layer", "dendrograms"] as const)
    expect(m[k], k).toEqual(want[k]);
  expect([...m.layers].sort((a, b) => a.name.localeCompare(b.name))).toEqual([...want.layers].sort((a, b) => a.name.localeCompare(b.name)));
  expect(Object.keys(m.obs).sort()).toEqual(Object.keys(want.obs).sort());
  expect(m.obs_kinds).toEqual(want.obs_kinds);
  for (const [k, col] of Object.entries(want.obs)) col.forEach((v, i) => close(m.obs[k][i], v, `${k}[${i}]`));
  for (const [k, s] of Object.entries(want._layer_sums)) {
    let t = 0;
    for (const v of got.layers[k]) if (!Number.isNaN(v)) t += v;
    expect(t, k).toBeCloseTo(s, 2);
  }
}

describe("browser readers", () => {
  it("read a cnvscope .h5ad like the server", async () => {
    same(await readH5ad(bytes("h5ad/demo_small.h5ad"), "demo_small"), expected("h5ad/demo_small.meta.json.gz"));
  });

  it("read an ASCAT.sc export with the same QC metrics as Python", async () => {
    same(await readAscatJson(bytes("ascatsc_json/demo_export.json.gz"), "demo_export"), expected("ascatsc_json/demo_export.meta.json.gz"));
  });
});
