import { useEffect, useState } from "react";
import { chrLabel } from "./canvas";

/** Bin index whose cumulative span holds gx (binary search). */
export function binAt(endCum: ArrayLike<number>, gx: number): number {
  let lo = 0, hi = endCum.length - 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (endCum[mid] < gx) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

/**
 * The bin drawn in each of `n` pixel columns across [gx0, gx1], or -1 for none.
 * A gap between two bins of the same chromosome (windows a tool filtered out,
 * e.g. epiAneufinder) takes the nearer bin, so the heatmap has no blank stripes;
 * gaps between chromosomes and beyond the last bin stay empty.
 */
export function binColumns(
  B: { chr: ArrayLike<string>; start_cum: ArrayLike<number>; end_cum: ArrayLike<number> },
  gx0: number, gx1: number, n: number, out: Float32Array | Int32Array = new Float32Array(n),
): Float32Array | Int32Array {
  for (let px = 0; px < n; px++) {
    const gx = gx0 + ((px + 0.5) / n) * (gx1 - gx0);
    const b = binAt(B.end_cum, gx);
    if (B.start_cum[b] <= gx && gx <= B.end_cum[b]) out[px] = b;
    else if (b > 0 && gx < B.start_cum[b] && B.chr[b - 1] === B.chr[b])
      out[px] = gx - B.end_cum[b - 1] < B.start_cum[b] - gx ? b - 1 : b;
    else out[px] = -1;
  }
  return out;
}

/** Visible bins [b0, b1) of the genomic window [gx0, gx1]. */
export function binRange(startCum: ArrayLike<number>, endCum: ArrayLike<number>, gx0: number, gx1: number): [number, number] {
  const n = endCum.length;
  let b0 = 0, b1 = n;
  for (let i = 0; i < n; i++) if (endCum[i] >= gx0) { b0 = i; break; }
  for (let i = n - 1; i >= 0; i--) if (startCum[i] <= gx1) { b1 = i + 1; break; }
  return [b0, b1];
}

/**
 * Chromosome labels that fit: [text, x] for spans [chr, left px, right px].
 * Each label needs its measured width plus a small gap from the previous one,
 * and the last chromosome keeps its label over the one before it.
 */
export function chromLabels(g: CanvasRenderingContext2D, spans: [string, number, number][]): [string, number][] {
  const out: [string, number, number][] = [], GAP = 2;
  for (const [chr, a, b] of spans) {
    const t = chrLabel(chr), w = g.measureText(t).width, x = (a + b) / 2;
    if (b - a < 4) continue;
    const prev = out[out.length - 1];
    if (prev && x - w / 2 < prev[1] + prev[2] / 2 + GAP) {
      if (chr !== spans[spans.length - 1][0]) continue;
      out.pop();
    }
    out.push([t, x, w]);
  }
  return out.map(([t, x]) => [t, x]);
}

/* ---------- cytobands, fetched lazily from the UCSC public API ---------- */
export interface Band { chr: string; start: number; end: number; stain: string }
export const STAIN_COLORS: Record<string, string> = {
  gneg: "#FFFFFF", gpos25: "#C8C8C8", gpos50: "#A0A0A0", gpos75: "#787878", gpos100: "#000000",
  acen: "#CC0000", gvar: "#A0A0A0", stalk: "#5C5C9C",
};
const cache = new Map<string, Band[]>();

export function useCytobands(genome: string | null): Band[] | null {
  const g = genome && /^(hg38|hg19|mm10|mm39)$/.test(genome) ? genome : null;
  const [bands, setBands] = useState<Band[] | null>(g ? cache.get(g) ?? null : null);
  useEffect(() => {
    if (!g) return;
    if (cache.has(g)) { setBands(cache.get(g)!); return; }
    let live = true;
    fetch(`https://api.genome.ucsc.edu/getData/track?genome=${g};track=cytoBandIdeo`)
      .then((r) => r.json())
      .then((json) => {
        type Raw = { chrom: string; chromStart: number; chromEnd: number; gieStain: string };
        // a list of bands, or (current API) bands grouped by chromosome
        const tr: Raw[] | Record<string, Raw[]> = json.cytoBandIdeo || [];
        const raw: Raw[] = Array.isArray(tr) ? tr : Object.values(tr).flat();
        const parsed = raw.map((b) => ({ chr: b.chrom, start: +b.chromStart, end: +b.chromEnd, stain: b.gieStain }));
        cache.set(g, parsed);
        if (live) setBands(parsed);
      })
      .catch(() => { /* offline: the ideogram is simply not drawn */ });
    return () => { live = false; };
  }, [g]);
  return bands;
}
