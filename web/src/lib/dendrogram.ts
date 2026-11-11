// Dendrogram layout from R hclust merge/height (ported from the ASCAT.sc
// Explorer). Leaves are -1-based row indices of the dataset; positive entries
// are 1-based earlier merges.

export interface Segment {
  x1: number;
  y1: number;
  x2: number;
  y2: number;
}

/**
 * Segments in (height, display row) space. `rowOf[leaf]` is the display row of
 * each dataset row, or -1 when the cell is hidden (its branch is pruned and the
 * parent keeps the surviving side connected).
 */
export function layoutDendrogram(merge: [number, number][], height: number[], rowOf: ArrayLike<number>): Segment[] {
  const segs: Segment[] = [];
  const nodeY = new Float64Array(merge.length).fill(NaN);
  const nodeX = new Float64Array(merge.length);
  const yOf = (idx: number): number | null => {
    if (idx < 0) {
      const r = rowOf[-idx - 1];
      return r != null && r >= 0 ? r + 0.5 : null;
    }
    const y = nodeY[idx - 1];
    return Number.isFinite(y) ? y : null;
  };
  const xOf = (idx: number) => (idx < 0 ? 0 : nodeX[idx - 1]);
  for (let i = 0; i < merge.length; i++) {
    const [l, r] = merge[i];
    const ly = yOf(l), ry = yOf(r);
    const h = height[i];
    nodeX[i] = h;
    if (ly == null && ry == null) continue;
    if (ly != null && ry != null) {
      nodeY[i] = (ly + ry) / 2;
      segs.push({ x1: xOf(l), y1: ly, x2: h, y2: ly });
      segs.push({ x1: xOf(r), y1: ry, x2: h, y2: ry });
      segs.push({ x1: h, y1: ly, x2: h, y2: ry });
    } else {
      const [c, y] = ly != null ? [l, ly] : [r, ry as number];
      nodeY[i] = y;
      segs.push({ x1: xOf(c), y1: y, x2: h, y2: y });
    }
  }
  return segs;
}

/**
 * Layout for the rows [r0, r1) on screen. When zoomed in on rows, only the
 * subtree joining the visible cells is drawn (their most recent common
 * ancestor), so its height scale fills the panel.
 */
export function layoutVisible(merge: [number, number][], height: number[], rowOf: ArrayLike<number>, r0: number, r1: number, nShown: number) {
  const full = () => ({ segments: layoutDendrogram(merge, height, rowOf), maxH: Math.max(0, ...height) });
  if (r0 <= 0 && r1 >= nShown) return full();
  const visible = new Set<number>();
  for (let i = 0; i < rowOf.length; i++) {
    const p = rowOf[i];
    if (p >= 0 && p >= r0 - 0.5 && p < r1 - 0.5) visible.add(i);
  }
  if (visible.size < 2) return full();
  const count = new Int32Array(merge.length);
  for (let i = 0; i < merge.length; i++) {
    const [l, r] = merge[i];
    count[i] = (l < 0 ? +visible.has(-l - 1) : count[l - 1]) + (r < 0 ? +visible.has(-r - 1) : count[r - 1]);
  }
  let root = merge.length - 1;
  for (let i = 0; i < merge.length; i++) if (count[i] === visible.size) { root = i; break; }
  const steps = new Set<number>();
  const stack = [root + 1];
  while (stack.length) {
    const code = stack.pop()!;
    if (code < 0) continue;
    steps.add(code - 1);
    const [l, r] = merge[code - 1];
    if (l > 0) stack.push(l);
    if (r > 0) stack.push(r);
  }
  const oldToNew = new Int32Array(merge.length).fill(-1);
  const m2: [number, number][] = [], h2: number[] = [];
  for (let i = 0; i < merge.length; i++) {
    if (!steps.has(i)) continue;
    oldToNew[i] = m2.length;
    const [l, r] = merge[i];
    m2.push([l < 0 ? l : oldToNew[l - 1] + 1, r < 0 ? r : oldToNew[r - 1] + 1]);
    h2.push(height[i]);
  }
  return { segments: layoutDendrogram(m2, h2, rowOf), maxH: Math.max(0, ...h2) };
}
