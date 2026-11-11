import { describe as group, expect, it } from "vitest";
import type { Meta } from "../api/client";
import { cellsTable, defaultFilters, describe, evaluate, fails, metricRanges } from "./filters";

const vals = {
  mapd: Float64Array.from([0.1, 0.2, 0.9, NaN]),
  ploidy: Float64Array.from([2, 2, 4, 1]),
  constant: Float64Array.from([1, 1, 1, 1]),
};

group("filters", () => {
  it("flags cells outside min and max, naming the reason", () => {
    const ev = evaluate(4, vals, { mapd: { min: null, max: 0.5 }, ploidy: { min: 1.5, max: 3 } });
    expect(ev.passing).toEqual([true, true, false, false]);
    expect(ev.reasons[2]).toBe("noisy, high ploidy");
    expect(ev.reasons[3]).toBe("low ploidy"); // NaN MAPD is not a failure
    expect(ev.failCounts).toEqual({ mapd: 1, ploidy: 2 });
    expect(ev.nPass).toBe(2);
  });
  it("ignores switched-off filters and metrics that do not vary", () => {
    expect(evaluate(4, vals, { mapd: { min: null, max: null } }).nPass).toBe(4);
    expect(metricRanges(vals).constant).toBeUndefined();
    expect(fails({ min: 1, max: null }, 0.5)).toBe(true);
    expect(fails({ min: 1, max: null }, NaN)).toBe(false);
  });
  it("sets the ASCAT.sc defaults only for metrics present", () => {
    const f = defaultFilters({ mapd: Float64Array.from([0.1, 0.2, 0.2, 0.21, 0.19, 0.9]) });
    expect(Object.keys(f)).toEqual(["mapd"]);
    expect(f.mapd.min).toBeNull();
    expect(f.mapd.max).toBeGreaterThan(0.2);
    // a default at the data's maximum flags nothing, so it is left out
    expect(defaultFilters({ bins_with_cna: Float64Array.from([1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3]) })).toEqual({});
  });
  it("describes and exports", () => {
    expect(describe("ploidy", { min: 1.5, max: 3 })).toBe("1.50 – 3.00");
    expect(describe("mapd", { min: null, max: 0.29 })).toBe("≤ 0.290");
    const meta = { version: "x", name: "d", tool: "t", n_cells: 4, cells: ["a", "b", "c", "d"], obs: { clone: ["A", "A", "B", null] }, obs_kinds: { clone: "categorical" } } as unknown as Meta;
    const f = { mapd: { min: null, max: 0.5 } };
    const tsv = cellsTable(meta, vals, f, evaluate(4, vals, f));
    const rows = tsv.trim().split("\n");
    expect(rows[2]).toBe("# filter: mapd ≤ 0.500");
    expect(rows[3].split("\t").slice(0, 4)).toEqual(["cell", "passes_qc", "failed", "clone"]);
    expect(rows[6].split("\t").slice(0, 4)).toEqual(["c", "FALSE", "noisy", "B"]);
    expect(rows[7].split("\t")[3]).toBe("");
  });
});
