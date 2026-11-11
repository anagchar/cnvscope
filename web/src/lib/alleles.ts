// Allele-specific states. ASCAT.sc's palette (get_allele_specific_colors() of ASCATsc_plot.R,
// hex values copied verbatim); the default palette and the state keys are in palettes.ts.

import type { Palette } from "./palettes";

export const AS_COLORS: Record<string, string> = {
  "0|0": "#153570", "1|0": "#345d8e", "1|1": "#c1c1c1", "2|0": "#737373", "2|1": "#edd78e",
  "3|0": "#c49631", "2|2": "#a8d9a8", "3|1": "#5fa85f", "4|0": "#2d7a2d", "3|2": "#f2b4a8",
  "4|1": "#cf3a3a", "5|0": "#7d1111", "3|3": "#d4b896", "4|2": "#a67c52", "5|1": "#7a5230",
  "6|0": "#4a2c12", "4|3": "#c7aed5", "5|2": "#9b6bb5", "6|1": "#6e3a8e", "7|0": "#4a1570",
  "7+": "#3d0d5c", NA: "#BEBEBE",
};

/** Colour index (into palette.asKeys) of every bin, NaN where missing or not a valid state. */
export function asStateMatrix(major: Float32Array, minor: Float32Array, palette: Palette): Float32Array {
  const index = new Map(palette.asKeys.map((k, i) => [k, i]));
  const out = new Float32Array(major.length);
  for (let i = 0; i < major.length; i++) {
    const k = palette.asKey(major[i], minor[i]);
    out[i] = k == null ? NaN : index.get(k)!;
  }
  return out;
}
