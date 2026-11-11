import { describe, expect, it } from "vitest";
import { madOutliers, madScaled, median, spearman } from "./qc";

describe("robust stats", () => {
  it("ignores non-finite values", () => {
    expect(median([3, NaN, 1, 2])).toBe(2);
    expect(madScaled([1, 2, 3, 4, 100])).toBeCloseTo(1.4826, 4);
  });

  it("flags only the bad side", () => {
    const x = [1, 1.1, 0.9, 1.05, 0.95, 10, -10];
    expect(madOutliers(x, "upper").flags).toEqual([false, false, false, false, false, true, false]);
    expect(madOutliers(x, "lower").flags).toEqual([false, false, false, false, false, false, true]);
    expect(madOutliers(x, "both").flags.filter(Boolean)).toHaveLength(2);
  });

  it("flags nothing when the MAD is zero", () => {
    expect(madOutliers([1, 1, 1, 5], "upper").flags.some(Boolean)).toBe(false);
  });
});

describe("spearman", () => {
  it("is 1 for monotone and -1 for reversed data", () => {
    expect(spearman([1, 2, 3, 4], [10, 20, 25, 100])).toBeCloseTo(1);
    expect(spearman([1, 2, 3, 4], [4, 3, 2, 1])).toBeCloseTo(-1);
  });
  it("needs three complete pairs", () => {
    expect(spearman([1, NaN, 3], [1, 2, 3])).toBeNaN();
  });
});
