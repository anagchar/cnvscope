// Robust statistics for QC thresholds (ported from the ASCAT.sc Explorer's qc.js).

export function quantile(sorted: ArrayLike<number>, p: number): number {
  if (!sorted.length) return NaN;
  const i = (sorted.length - 1) * p, a = Math.floor(i), b = Math.ceil(i);
  return sorted[a] + (sorted[b] - sorted[a]) * (i - a);
}

const finiteSorted = (x: ArrayLike<number>) => Array.from(x).filter((v) => Number.isFinite(v)).sort((a, b) => a - b);

export function median(x: ArrayLike<number>): number {
  return quantile(finiteSorted(x), 0.5);
}

/** MAD scaled to the SD of a normal distribution. */
export function madScaled(x: ArrayLike<number>): number {
  const v = finiteSorted(x);
  if (!v.length) return NaN;
  const m = quantile(v, 0.5);
  return 1.4826 * quantile(v.map((a) => Math.abs(a - m)).sort((a, b) => a - b), 0.5);
}

export type BadSide = "upper" | "lower" | "both";

/** Thresholds at median ± k·MAD and a per-value flag on the metric's bad side. */
export function madOutliers(x: ArrayLike<number>, dir: BadSide, k = 3) {
  const med = median(x), mad = madScaled(x);
  const lo = med - k * mad, hi = med + k * mad;
  const flags = Array.from(x, (v) =>
    !Number.isFinite(v) || !(mad > 0) ? false : (dir !== "lower" && v > hi) || (dir !== "upper" && v < lo),
  );
  return { median: med, mad, lo, hi, flags };
}

function corr(a: number[], b: number[]): number {
  const n = a.length;
  let ma = 0, mb = 0;
  for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < n; i++) {
    const da = a[i] - ma, db = b[i] - mb;
    sab += da * db; saa += da * da; sbb += db * db;
  }
  return saa > 0 && sbb > 0 ? sab / Math.sqrt(saa * sbb) : NaN;
}

/** Spearman rank correlation over pairs where both values are finite; ties get mean ranks. */
export function spearman(a: ArrayLike<number>, b: ArrayLike<number>): number {
  const idx: number[] = [];
  for (let i = 0; i < a.length; i++) if (Number.isFinite(a[i]) && Number.isFinite(b[i])) idx.push(i);
  if (idx.length < 3) return NaN;
  const rank = (v: ArrayLike<number>) => {
    const o = idx.map((i, j) => [v[i], j] as const).sort((p, q) => p[0] - q[0]);
    const r = new Array<number>(o.length);
    for (let s = 0; s < o.length; ) {
      let e = s;
      while (e + 1 < o.length && o[e + 1][0] === o[s][0]) e++;
      for (let k = s; k <= e; k++) r[o[k][1]] = (s + e) / 2;
      s = e + 1;
    }
    return r;
  };
  return corr(rank(a), rank(b));
}
