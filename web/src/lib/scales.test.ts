import { describe, expect, it } from "vitest";
import { buildLut, LINEAR_LUT, scaleFor, valueToUnit } from "./scales";

describe("scaleFor", () => {
  it("maps integer copy number to discrete levels with 6+ clamped", () => {
    const s = scaleFor("integer_cn");
    expect(s.kind).toBe("discrete");
    expect(s.baseline).toBe(2);
    expect(valueToUnit(s, 0)).toBe(0);
    expect(valueToUnit(s, 2.4)).toBe(2);
    expect(s.labels.at(-1)).toBe("11+");
    expect(valueToUnit(s, 11)).toBe(11);
    expect(valueToUnit(s, 40)).toBe(11);
    expect(valueToUnit(s, -3)).toBe(0);
    expect(valueToUnit(s, NaN)).toBeNaN();
  });

  it("maps ternary states to loss/neutral/gain", () => {
    const s = scaleFor("ternary_state");
    expect([-1, 0, 1].map((v) => valueToUnit(s, v))).toEqual([0, 1, 2]);
    expect(s.labels).toEqual(["Loss", "Normal", "Gain"]);
  });

  it("centres log ratios on 0 and clamps", () => {
    const s = scaleFor("continuous_ratio");
    expect(valueToUnit(s, 0)).toBe(0.5);
    expect(valueToUnit(s, -5)).toBe(0);
    expect(valueToUnit(s, 5)).toBe(1);
  });

  it("uses the data maximum for counts", () => {
    expect(scaleFor("count", { countMax: 40 }).hi).toBe(40);
    expect(scaleFor("count").hi).toBe(1);
  });
});

describe("buildLut", () => {
  const resolve = (t: string): [number, number, number] =>
    t === "--a" ? [0, 0, 0] : [255, 255, 255];

  it("has one entry per discrete level", () => {
    const lut = buildLut({ kind: "discrete", lo: 0, hi: 1, stops: ["--a", "--b"], labels: [], baseline: 0 }, resolve);
    expect(Array.from(lut)).toEqual([0, 0, 0, 255, 255, 255, 255, 255]);
  });

  it("interpolates linear stops end to end", () => {
    const lut = buildLut({ kind: "linear", lo: 0, hi: 1, stops: ["--a", "--b"], labels: [], baseline: 0 }, resolve);
    expect(lut.length).toBe(LINEAR_LUT * 4);
    expect(Array.from(lut.slice(0, 4))).toEqual([0, 0, 0, 255]);
    expect(Array.from(lut.slice(-4))).toEqual([255, 255, 255, 255]);
  });
});
