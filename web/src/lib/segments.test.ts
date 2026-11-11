import { describe, expect, it } from "vitest";
import { asStateMatrix } from "./alleles";
import { PALETTES } from "./palettes";
import { binsAltered, filterRow, filterSmallSegments, followFilter } from "./segments";

describe("filterSmallSegments", () => {
  it("resets short non-baseline runs to the baseline", () => {
    const row = Float32Array.from([2, 3, 3, 2, 1, 1, 1, 2, 4]);
    filterRow(row, 3, 2);
    expect(Array.from(row)).toEqual([2, 2, 2, 2, 1, 1, 1, 2, 2]);
  });
  it("leaves missing bins alone and works per row", () => {
    const m = Float32Array.from([3, NaN, 3, 3, 3, 3]);
    const f = filterSmallSegments(m, 2, 3, 2, 2);
    expect(Array.from(f)).toEqual([2, NaN, 2, 3, 3, 3]);
    expect(Array.from(m)).toEqual([3, NaN, 3, 3, 3, 3]); // input untouched
  });
  it("moves the alleles back to 1|1 where the total was reset", () => {
    const total = Float32Array.from([3, 2]), filtered = Float32Array.from([2, 2]);
    expect(Array.from(followFilter(Float32Array.from([2, 1]), total, filtered, 2, 1))).toEqual([1, 1]);
  });
  it("counts altered bins", () => {
    expect(Array.from(binsAltered(Float32Array.from([2, 3, NaN, 1, 2, 2]), 2, 3, 2))).toEqual([1, 1]);
  });
});

describe("allele states", () => {
  it("orders major|minor and lumps high states", () => {
    const { ascat, signals } = PALETTES;
    expect(ascat.asKey(1, 2)).toBe("2|1");
    expect(ascat.asKey(5, 3)).toBe("7+");
    expect(ascat.asKey(NaN, 1)).toBeNull();
    expect(ascat.asKey(2, -1)).toBeNull(); // not a state: NA, never 7+
    // SIGNALS: states up to 4 copies, then by total
    expect(signals.asKey(2, 0)).toBe("2|0");
    expect(signals.asKey(3, 2)).toBe("5");
    expect(signals.asKey(9, 4)).toBe("11+");
    expect(signals.asKey(2, -1)).toBeNull();
    const m = asStateMatrix(Float32Array.from([1, 2]), Float32Array.from([1, NaN]), signals);
    expect(signals.asKeys[m[0]]).toBe("1|1");
    expect(m[1]).toBeNaN();
  });
});
