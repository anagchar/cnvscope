// Leaf order of an R hclust merge matrix (left child first), as row indices.
// A port of src/cnvscope/analysis.py: hclust_order.
export function layoutOrder(merge: [number, number][]): number[] {
  if (!merge.length) return [0];
  const out: number[] = [], stack = [merge.length];
  while (stack.length) {
    const code = stack.pop()!;
    if (code < 0) { out.push(-code - 1); continue; }
    const [left, right] = merge[code - 1];
    stack.push(right, left);
  }
  return out;
}
