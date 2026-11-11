// Colour palettes, taken from the single-cell copy-number literature.
//
// "signals" (default): the palette of the Shah lab's DLP+ / SIGNALS / scgenome tools
//   (Laks et al. 2019 Cell; Funnell et al. 2022 Nature; shahcompbio/signals R/col_palettes.R),
//   built on ColorBrewer ramps: blues for losses, grey for 2 copies, orange → red → magenta for
//   gains. Allele-specific states follow scCNAS_colors (copy-neutral LOH 2|0 in dark grey,
//   states above 4 copies by their total), BAF the A-Hom / balanced / B-Hom ramp (green–grey–
//   purple). It is the de facto standard of single-cell CNV figures and safe for the common
//   colour-vision deficiencies.
// "ascat": ASCAT.sc's own colours (get_allele_specific_colors() of ASCATsc_plot.R), to match
//   plots made in R.
//
// Continuous scales use evenly-graded maps (Crameri, Shephard & Heron 2020, Nat Commun 11:5444):
// a blue–grey–red diverging map for ratios (matching the copy-number colours) and viridis for
// counts and probabilities; categorical annotations use the Okabe–Ito colour-blind-safe set.

import { AS_COLORS } from "./alleles";

export type PaletteName = "signals" | "ascat";

export interface Palette {
  name: PaletteName;
  label: string;
  /** colours of total copy number 0 … n-1, the last meaning "n-1 or more" */
  cn: string[];
  /** allele-specific states in legend order, and their colours */
  asKeys: string[];
  asColors: string[];
  /** state key of a major/minor pair, null when missing or invalid */
  asKey: (major: number, minor: number) => string | null;
}

const SIGNALS_CN = ["#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"];
const SIGNALS_AS: [string, string][] = [
  ["0|0", "#3182BD"], ["1|0", "#9ECAE1"], ["1|1", "#CCCCCC"], ["2|0", "#666666"], ["2|1", "#FDCC8A"],
  ["3|0", "#FEE2BC"], ["2|2", "#FC8D59"], ["3|1", "#FDC1A4"], ["4|0", "#FB590E"],
  ["5", "#E34A33"], ["6", "#B30000"], ["7", "#980043"], ["8", "#DD1C77"], ["9", "#DF65B0"], ["10", "#C994C7"], ["11+", "#D4B9DA"],
];
const SIGNALS_AS_SET = new Set(SIGNALS_AS.map(([k]) => k));

function ordered(a: number, b: number): [number, number] | null {
  if (!Number.isFinite(a) || !Number.isFinite(b)) return null;
  const M = Math.round(Math.max(a, b)), m = Math.round(Math.min(a, b));
  return m < 0 ? null : [M, m];
}

export const PALETTES: Record<PaletteName, Palette> = {
  signals: {
    name: "signals",
    label: "SIGNALS",
    cn: SIGNALS_CN,
    asKeys: SIGNALS_AS.map(([k]) => k),
    asColors: SIGNALS_AS.map(([, c]) => c),
    asKey: (a, b) => {
      const p = ordered(a, b);
      if (!p) return null;
      const [M, m] = p, t = M + m;
      if (t >= 11) return "11+";
      if (t >= 5) return String(t);
      const k = `${M}|${m}`;
      return SIGNALS_AS_SET.has(k) ? k : null;
    },
  },
  ascat: (() => {
    const keys = Object.keys(AS_COLORS).filter((k) => k !== "NA");
    const set = new Set(keys);
    return {
      name: "ascat" as const,
      label: "ASCAT.sc",
      cn: ["#1c5cab", "#86b6ef", "#f0efec", "#f4a59c", "#e8604f", "#a82628", "#5a0c10"],
      asKeys: keys,
      asColors: keys.map((k) => AS_COLORS[k]),
      asKey: (a: number, b: number) => {
        const p = ordered(a, b);
        if (!p) return null;
        const [M, m] = p;
        if (M + m > 7) return "7+";
        const k = `${M}|${m}`;
        return set.has(k) ? k : null;
      },
    };
  })(),
};

/** Ratios: blue (loss) – grey – red (gain), the copy-number colours' ends. */
export const DIVERGING = ["#3182BD", "#9ECAE1", "#CCCCCC", "#FC8D59", "#B30000"];
/** viridis (van der Walt & Smith), evenly graded and colour-blind safe. */
export const SEQUENTIAL = ["#440154", "#3B528B", "#21918C", "#5EC962", "#FDE725"];
/** BAF: A-Hom (green) – balanced (grey) – B-Hom (purple), from SIGNALS. */
export const BAF = ["#56941E", "#94C773", "#D5D5D4", "#7B52AE", "#471871"];
/** Okabe & Ito (2008), then Paul Tol's "muted" set for more categories. */
export const CATEGORICAL = [
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999",
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499",
];

/**
 * On a dark background the neutral state is darkened so most of the heatmap doesn't glare,
 * and copy-neutral LOH (dark grey in SIGNALS) turns light grey so it stays apart from it.
 */
export function forTheme(colors: string[], dark: boolean, neutral: string): string[] {
  if (!dark) return colors;
  return colors.map((c) => {
    const u = c.toUpperCase();
    if (u === neutral.toUpperCase()) return "#46464C";
    if (u === "#666666") return "#A9A9B1";
    return c;
  });
}
