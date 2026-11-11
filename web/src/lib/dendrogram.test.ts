import { describe, expect, it } from "vitest";
import { layoutDendrogram } from "./dendrogram";

// leaves 0,1 merge at 1; leaf 2 joins at 2
const merge: [number, number][] = [[-1, -2], [-3, 1]];
const height = [1, 2];

describe("layoutDendrogram", () => {
  it("draws three segments per merge, centred between children", () => {
    const segs = layoutDendrogram(merge, height, [0, 1, 2]);
    expect(segs).toHaveLength(6);
    // vertical of the first merge spans rows 0.5..1.5 at height 1
    expect(segs[2]).toEqual({ x1: 1, y1: 0.5, x2: 1, y2: 1.5 });
    // second merge joins leaf 2 (y 2.5) with the first merge's midpoint (y 1)
    expect(segs[5]).toEqual({ x1: 2, y1: 2.5, x2: 2, y2: 1 });
  });

  it("keeps the surviving branch when a leaf is hidden", () => {
    const segs = layoutDendrogram(merge, height, [0, -1, 1]);
    // merge 0 has one visible child: one horizontal segment
    expect(segs[0]).toEqual({ x1: 0, y1: 0.5, x2: 1, y2: 0.5 });
    expect(segs).toHaveLength(4);
  });
});
