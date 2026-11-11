// Colour scales chosen by what a layer's values mean, never by the tool. Colours come
// from the active palette (palettes.ts); stops are literal colours or theme tokens.

import type { Semantics } from "../api/client";
import { BAF, DIVERGING, PALETTES, SEQUENTIAL, forTheme, type Palette } from "./palettes";

export interface Scale {
  /** discrete: one colour per integer level lo..hi (values clamped); linear: stops spread over lo..hi */
  kind: "discrete" | "linear";
  lo: number;
  hi: number;
  /** colours (or CSS custom property names) */
  stops: string[];
  /** level labels for the legend (discrete) or end labels (linear) */
  labels: string[];
  /** the value that means "no change" */
  baseline: number;
}

export interface ScaleOpts {
  countMax?: number;
  palette?: Palette;
  dark?: boolean;
}

export function scaleFor(sem: Semantics, opts: ScaleOpts = {}): Scale {
  const pal = opts.palette ?? PALETTES.signals;
  const dark = !!opts.dark;
  const cn = forTheme(pal.cn, dark, pal.cn[2]);
  switch (sem) {
    case "integer_cn": {
      const n = cn.length;
      const labels = cn.map((_, i) => (i === n - 1 ? `${i}+` : String(i)));
      return { kind: "discrete", lo: 0, hi: n - 1, stops: cn, labels, baseline: 2 };
    }
    case "allele_cn":
      return { kind: "discrete", lo: 0, hi: 4, stops: [cn[0], cn[2], cn[3], cn[4], cn[5]], labels: ["0", "1", "2", "3", "4+"], baseline: 1 };
    case "ternary_state":
      return { kind: "discrete", lo: -1, hi: 1, stops: [cn[0], cn[2], cn[5]], labels: ["Loss", "Normal", "Gain"], baseline: 0 };
    case "continuous_ratio":
      return { kind: "linear", lo: -1, hi: 1, stops: forTheme(DIVERGING, dark, DIVERGING[2]), labels: ["≤ −1", "0", "≥ 1"], baseline: 0 };
    case "continuous_cn":
      return { kind: "linear", lo: 0, hi: 6, stops: cn.slice(0, 7), labels: ["0", "2", "6+"], baseline: 2 };
    case "ratio":
      return { kind: "linear", lo: 0, hi: 2, stops: forTheme(DIVERGING, dark, DIVERGING[2]), labels: ["0", "1", "≥ 2"], baseline: 1 };
    case "baf":
      return { kind: "linear", lo: 0, hi: 1, stops: forTheme(BAF, dark, BAF[2]), labels: ["0", "0.5", "1"], baseline: 0.5 };
    case "probability":
      return { kind: "linear", lo: 0, hi: 1, stops: SEQUENTIAL, labels: ["0", "1"], baseline: 0 };
    case "count": {
      const hi = opts.countMax && opts.countMax > 0 ? opts.countMax : 1;
      return { kind: "linear", lo: 0, hi, stops: SEQUENTIAL, labels: ["0", String(Math.round(hi))], baseline: 0 };
    }
  }
}

/** Major|minor states of the palette (values are indices into palette.asKeys). */
export function alleleScale(palette: Palette = PALETTES.signals, dark = false): Scale {
  const neutral = palette.asColors[palette.asKeys.indexOf("1|1")];
  return {
    kind: "discrete",
    lo: 0,
    hi: palette.asKeys.length - 1,
    stops: forTheme(palette.asColors, dark, neutral),
    labels: palette.asKeys,
    baseline: palette.asKeys.indexOf("1|1"),
  };
}

/** Colour of a value on a discrete scale. */
export function discreteColor(s: Scale, v: number): string | null {
  if (!Number.isFinite(v)) return null;
  return s.stops[Math.min(s.hi, Math.max(s.lo, Math.round(v))) - s.lo];
}

/** Position of a value on the scale: the colour index (discrete) or 0..1 (linear); NaN stays NaN. */
export function valueToUnit(s: Scale, v: number): number {
  if (Number.isNaN(v)) return NaN;
  if (s.kind === "discrete") return Math.min(s.hi, Math.max(s.lo, Math.round(v))) - s.lo;
  return Math.min(1, Math.max(0, (v - s.lo) / (s.hi - s.lo)));
}

export const LINEAR_LUT = 256;

/** RGBA lookup table: one entry per level (discrete) or LINEAR_LUT entries (linear). */
export function buildLut(s: Scale, resolve: (token: string) => [number, number, number]): Uint8Array {
  const cols = s.stops.map(resolve);
  if (s.kind === "discrete") {
    const out = new Uint8Array(cols.length * 4);
    cols.forEach((c, i) => out.set([c[0], c[1], c[2], 255], i * 4));
    return out;
  }
  const out = new Uint8Array(LINEAR_LUT * 4);
  for (let i = 0; i < LINEAR_LUT; i++) {
    const t = (i / (LINEAR_LUT - 1)) * (cols.length - 1);
    const a = Math.min(cols.length - 2, Math.floor(t));
    const f = t - a;
    const [c0, c1] = [cols[a], cols[a + 1]];
    out.set([0, 1, 2].map((k) => Math.round(c0[k] + (c1[k] - c0[k]) * f)).concat(255), i * 4);
  }
  return out;
}

/** Whether a layer holds calls (drawn as steps) rather than measurements (drawn as points). */
export const isCall = (sem: Semantics) => sem === "integer_cn" || sem === "allele_cn" || sem === "ternary_state";
