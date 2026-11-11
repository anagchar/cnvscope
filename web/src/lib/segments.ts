// Minimum segment size: runs of a non-baseline call shorter than minBins are shown
// as the baseline (ported from the ASCAT.sc Explorer's filterSmallSegments).

/** Filter one row in place; returns whether anything changed. */
export function filterRow(row: Float32Array, minBins: number, baseline: number): boolean {
  let changed = false;
  for (let i = 0; i < row.length; ) {
    if (!Number.isFinite(row[i])) { i++; continue; }
    const cn = Math.round(row[i]);
    if (cn === baseline) { i++; continue; }
    let j = i + 1;
    while (j < row.length && Number.isFinite(row[j]) && Math.round(row[j]) === cn) j++;
    if (j - i < minBins) {
      row.fill(baseline, i, j);
      changed = true;
    }
    i = j;
  }
  return changed;
}

/** A filtered copy of a rows × cols matrix. */
export function filterSmallSegments(m: Float32Array, rows: number, cols: number, minBins: number, baseline: number): Float32Array {
  const out = m.slice();
  if (minBins <= 1) return out;
  for (let r = 0; r < rows; r++) filterRow(out.subarray(r * cols, (r + 1) * cols), minBins, baseline);
  return out;
}

/**
 * Allele layers following a filtered total: where the total was reset to the
 * baseline, the alleles become `alleleBaseline` (1|1 for a diploid baseline).
 */
export function followFilter(allele: Float32Array, total: Float32Array, filtered: Float32Array, baseline: number, alleleBaseline: number): Float32Array {
  const out = allele.slice();
  for (let i = 0; i < out.length; i++)
    if (Number.isFinite(total[i]) && Math.round(total[i]) !== baseline && Math.round(filtered[i]) === baseline) out[i] = alleleBaseline;
  return out;
}

/** Bins per row whose call is not the baseline. */
export function binsAltered(m: Float32Array, rows: number, cols: number, baseline: number): Float64Array {
  const out = new Float64Array(rows);
  for (let r = 0; r < rows; r++) {
    let n = 0;
    for (let c = r * cols, e = c + cols; c < e; c++) if (Number.isFinite(m[c]) && Math.round(m[c]) !== baseline) n++;
    out[r] = n;
  }
  return out;
}
