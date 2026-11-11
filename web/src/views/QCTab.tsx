import { useEffect, useMemo, useRef, useState, type ReactNode } from "react";
import { FONT, axes, axisTicks, css, emptyMsg, esc, fit, legendRow, mixRGB, niceTicks, rgbOf, roundTop, tip } from "../lib/canvas";
import { cellsTable, defaultFilters, describe, downloadText, fails, isActive, type Filter } from "../lib/filters";
import { INT_METRICS, fmtMetric, type Metric } from "../lib/metrics";
import { DIVERGING } from "../lib/palettes";
import { spearman } from "../lib/qc";
import { useDerived } from "../state/derived";
import { useStore } from "../state/store";
import { SegmentControl } from "./SegmentControl";
import { DlIcon, useSize, useWidth } from "./widgets";

type Vals = Record<string, Float64Array>;

function statusHTML(pass: boolean, reason: string) {
  return pass
    ? `<span class="status" style="color:var(--good)">● passes QC</span>`
    : `<span class="status" style="color:var(--serious)">▲ flagged${reason ? ": " + esc(reason) : ""}</span>`;
}

/* ---------- Explore: any metric against any other ---------- */
type Range = [number, number];
interface View { x: Range; y: Range }

/**
 * Drag a box to zoom, scroll to zoom at the cursor, shift-drag to pan, double-click to reset.
 * The QC filter lines can be dragged: they set the filters (same as on the histograms).
 */
function Scatter({ cells, mx, my, vals, passing, shown, reasons, sel, onSelect }: {
  cells: string[]; mx?: Metric; my?: Metric; vals: Vals; passing: boolean[]; shown: Set<number>; reasons: string[];
  sel: number | null; onSelect: (i: number) => void;
}) {
  useStore((s) => s.theme); // redraw (every render) when the theme changes
  const { ranges, filters } = useDerived();
  const setFilter = useStore((s) => s.setFilter);
  const wrap = useRef<HTMLDivElement>(null), cv = useRef<HTMLCanvasElement>(null);
  const pts = useRef<[number, number, number][]>([]), hover = useRef<number | null>(null);
  const geo = useRef<{ tx: (v: number) => number; ty: (v: number) => number; ix: (p: number) => number; iy: (p: number) => number; box: { l: number; t: number; w: number; h: number } } | null>(null);
  const { w: W, h: H } = useSize(wrap, 340);
  const [view, setView] = useState<View | null>(null);
  const drag = useRef<
    | { kind: "box" | "pan"; sx: number; sy: number; cx: number; cy: number; v0: View }
    | { kind: "thresh"; metric: string; side: "min" | "max"; isX: boolean; lim: Range }
    | null
  >(null);
  const [, redraw] = useState(0);
  useEffect(() => setView(null), [mx?.key, my?.key]);

  // the threshold (sidebar key) drawn for a metric, if one is active
  // the active filter bounds of a metric: [side, value]
  const boundsOf = (m?: Metric): ["min" | "max", number][] => {
    const f = m && filters[m.key];
    if (!f || !isActive(f) || !ranges[m.key]) return [];
    return ([["min", f.min], ["max", f.max]] as const).filter((b): b is ["min" | "max", number] => b[1] != null);
  };

  // full data range of an axis, padded
  const auto = useMemo(() => {
    if (!mx || !my) return null;
    const range = (v: Float64Array, log: boolean): Range | null => {
      const ok = Array.from(v).filter((x) => !Number.isNaN(x) && (!log || x > 0));
      if (!ok.length) return null;
      let a = Math.min(...ok), b = Math.max(...ok);
      if (a === b) { a -= log ? a * 0.5 : 0.5; b += log ? b * 0.5 : 0.5; }
      if (log) return [a / 1.15, b * 1.15];
      const pad = (b - a) * 0.06;
      return [a - pad, b + pad];
    };
    const x = range(vals[mx.key], mx.log), y = range(vals[my.key], my.log);
    return x && y ? { x, y } : null;
  }, [mx, my, vals]);

  useEffect(() => {
    const g = fit(cv.current!, W, H);
    g.clearRect(0, 0, W, H);
    geo.current = null;
    if (!mx || !my) return emptyMsg(g, W, H, "No metric varies between the cells.");
    if (!auto) return emptyMsg(g, W, H, "No cell has both values.");
    const xs = vals[mx.key], ys = vals[my.key];
    const [x0, x1] = view?.x ?? auto.x, [y0, y1] = view?.y ?? auto.y;
    const box = { l: 52, t: 30, w: W - 66, h: H - 74 };
    const tk = (a: number, b: number, log: boolean) => (log ? axisTicks(a, b) : niceTicks(a, b, 5));
    const f = (m: Metric) => (v: number) => fmtMetric(m.key, v).replace(/\.0+(?=[kM]?$)|(\.\d*?)0+(?=[kM]?$)/, "$1");
    const { tx, ty } = axes(g, box, tk(x0, x1, mx.log), tk(y0, y1, my.log), { x0, x1, y0, y1, xlog: mx.log, ylog: my.log, xlab: mx.label, ylab: my.label, xfmt: f(mx), yfmt: f(my) });
    const inv = (p: number, a: number, b: number, log: boolean) => (log ? Math.pow(10, Math.log10(a) + p * (Math.log10(b) - Math.log10(a))) : a + p * (b - a));
    geo.current = { tx, ty, box, ix: (px) => inv((px - box.l) / box.w, x0, x1, mx.log), iy: (py) => inv((box.t + box.h - py) / box.h, y0, y1, my.log) };
    g.save();
    g.beginPath(); g.rect(box.l, box.t, box.w, box.h); g.clip();
    g.setLineDash([4, 4]); g.strokeStyle = css("--ink3"); g.lineWidth = 1;
    const line = (v: number, isX: boolean) => {
      const lim = isX ? [x0, x1] : [y0, y1];
      if (!(Number.isFinite(v) && v > lim[0] && v < lim[1])) return false;
      g.beginPath();
      if (isX) { g.moveTo(tx(v), box.t); g.lineTo(tx(v), box.t + box.h); } else { g.moveTo(box.l, ty(v)); g.lineTo(box.l + box.w, ty(v)); }
      g.stroke();
      return true;
    };
    for (const [m, isX] of [[mx, true], [my, false]] as const)
      for (const [v, side] of [[m.lo, "lower"], [m.hi, "upper"]] as const) if (m.dir === side || m.dir === "both") line(v, isX);
    g.setLineDash([]); g.strokeStyle = css("--qc-fail"); g.globalAlpha = 0.85;
    let drawnThreshold = false;
    for (const [m, isX] of [[mx, true], [my, false]] as const)
      for (const [side, v] of boundsOf(m)) {
        const dr = drag.current;
        g.lineWidth = dr?.kind === "thresh" && dr.isX === isX && dr.side === side ? 2.5 : 1.5;
        if (line(v, isX)) drawnThreshold = true;
      }
    g.globalAlpha = 1;
    pts.current = [];
    const inside = (x: number, y: number) => x >= box.l - 4 && x <= box.l + box.w + 4 && y >= box.t - 4 && y <= box.t + box.h + 4;
    for (const pass of [true, false])
      for (let i = 0; i < cells.length; i++) {
        if (passing[i] !== pass || Number.isNaN(xs[i]) || Number.isNaN(ys[i]) || (mx.log && !(xs[i] > 0)) || (my.log && !(ys[i] > 0))) continue;
        const x = tx(xs[i]), y = ty(ys[i]);
        if (!inside(x, y)) continue;
        pts.current.push([x, y, i]);
        g.globalAlpha = shown.has(i) || !pass ? 0.85 : 0.18;
        g.fillStyle = css(pass ? "--qc-pass" : "--qc-fail");
        g.beginPath();
        if (pass) g.arc(x, y, 3.2, 0, 7);
        else { g.moveTo(x, y - 4); g.lineTo(x + 3.6, y + 2.6); g.lineTo(x - 3.6, y + 2.6); }
        g.fill();
      }
    g.globalAlpha = 1;
    if (sel != null && !Number.isNaN(xs[sel]) && !Number.isNaN(ys[sel])) {
      const x = tx(xs[sel]), y = ty(ys[sel]);
      g.strokeStyle = css("--panel"); g.lineWidth = 2; g.beginPath(); g.arc(x, y, 7, 0, 7); g.stroke();
      g.strokeStyle = css("--accent"); g.beginPath(); g.arc(x, y, 5.5, 0, 7); g.stroke();
    }
    const d = drag.current;
    if (d?.kind === "box") {
      g.fillStyle = css("--accent-soft"); g.strokeStyle = css("--accent"); g.lineWidth = 1;
      const bx = Math.min(d.sx, d.cx), by = Math.min(d.sy, d.cy);
      g.fillRect(bx, by, Math.abs(d.cx - d.sx), Math.abs(d.cy - d.sy));
      g.strokeRect(bx, by, Math.abs(d.cx - d.sx), Math.abs(d.cy - d.sy));
    }
    g.restore();
    const items: [string, string, "dot" | "tri" | "line"][] = [["--qc-pass", "passes QC", "dot"], ["--qc-fail", "flagged", "tri"], ["--ink3", "median ± 3 MAD", "line"]];
    if (drawnThreshold) items.push(["--qc-fail", "filter (drag it)", "line"]);
    legendRow(g, box.l, 12, items);
  });

  const at = (e: { clientX: number; clientY: number }) => {
    const r = cv.current!.getBoundingClientRect();
    return { x: e.clientX - r.left, y: e.clientY - r.top };
  };
  // a threshold line within a few pixels of the pointer
  const lineAt = (x: number, y: number) => {
    const G = geo.current;
    if (!G || !mx || !my) return null;
    const { box } = G;
    if (x < box.l || x > box.l + box.w || y < box.t || y > box.t + box.h) return null;
    for (const [m, isX] of [[mx, true], [my, false]] as const)
      for (const [side, v] of boundsOf(m)) {
        const p = isX ? G.tx(v) : G.ty(v);
        if (Math.abs((isX ? x : y) - p) <= 5) return { metric: m.key, side, isX, lim: ranges[m.key] };
      }
    return null;
  };
  const cur = () => (view ?? auto)!;
  const onDown = (e: React.MouseEvent) => {
    if (e.button !== 0 || !auto) return;
    const p = at(e);
    const l = lineAt(p.x, p.y);
    tip(null);
    if (l) drag.current = { kind: "thresh", ...l };
    else drag.current = { kind: e.shiftKey ? "pan" : "box", sx: p.x, sy: p.y, cx: p.x, cy: p.y, v0: cur() };
  };
  const onMove = (e: React.MouseEvent) => {
    const p = at(e), d = drag.current, G = geo.current;
    if (d && G) {
      if (d.kind === "thresh") {
        let v = d.isX ? G.ix(p.x) : G.iy(p.y);
        v = Math.max(d.lim[0], Math.min(d.lim[1], v));
        v = INT_METRICS.has(d.metric) ? Math.round(v) : +v.toPrecision(4);
        const f: Filter = useStore.getState().filters?.[d.metric] ?? { min: null, max: null };
        setFilter(d.metric, { ...f, [d.side]: v });
      } else if (d.kind === "pan") {
        const shift = (r: Range, dp: number, len: number, log: boolean): Range => {
          if (log) { const k = Math.pow(r[1] / r[0], -dp / len); return [r[0] * k, r[1] * k]; }
          const s = (-dp / len) * (r[1] - r[0]);
          return [r[0] + s, r[1] + s];
        };
        setView({ x: shift(d.v0.x, p.x - d.sx, G.box.w, !!mx?.log), y: shift(d.v0.y, d.sy - p.y, G.box.h, !!my?.log) });
      } else {
        d.cx = p.x; d.cy = p.y;
        redraw((n) => n + 1);
      }
      return;
    }
    const l = lineAt(p.x, p.y);
    let best: number | null = null, bd = 64;
    if (!l) for (const [px, py, i] of pts.current) { const dd = (px - p.x) ** 2 + (py - p.y) ** 2; if (dd < bd) { bd = dd; best = i; } }
    cv.current!.style.cursor = l ? (l.isX ? "ew-resize" : "ns-resize") : best == null ? "crosshair" : "pointer";
    hover.current = best;
    tip(e, l ? `drag to move the <b>${l.side === "max" ? "maximum" : "minimum"} ${esc((l.isX ? mx : my)!.label.toLowerCase())}</b> filter`
      : best == null || !mx || !my ? null
      : `<b class="mono">${esc(cells[best])}</b><br>${statusHTML(passing[best], reasons[best])}<br><span class="k">${esc(mx.label)}</span> ${fmtMetric(mx.key, vals[mx.key][best])}<br><span class="k">${esc(my.label)}</span> ${fmtMetric(my.key, vals[my.key][best])}`);
  };
  const onUp = (e: React.MouseEvent) => {
    const d = drag.current, G = geo.current;
    drag.current = null;
    if (!d || !G) return;
    if (d.kind === "box") {
      const p = at(e);
      if (Math.abs(p.x - d.sx) < 5 && Math.abs(p.y - d.sy) < 5) {
        if (hover.current != null) onSelect(hover.current);
      } else if (Math.abs(p.x - d.sx) > 8 && Math.abs(p.y - d.sy) > 8) {
        const xa = G.ix(Math.min(d.sx, p.x)), xb = G.ix(Math.max(d.sx, p.x));
        const ya = G.iy(Math.max(d.sy, p.y)), yb = G.iy(Math.min(d.sy, p.y));
        setView({ x: [xa, xb], y: [ya, yb] });
      }
    }
    redraw((n) => n + 1);
  };
  useEffect(() => {
    const el = cv.current;
    if (!el) return;
    const h = (e: WheelEvent) => {
      const G = geo.current, v = view ?? auto;
      if (!G || !v || !mx || !my) return;
      const p = at(e);
      if (p.x < G.box.l || p.x > G.box.l + G.box.w || p.y < G.box.t || p.y > G.box.t + G.box.h) return;
      e.preventDefault();
      const k = Math.exp(e.deltaY * 0.0015);
      const zoom = (r: Range, c: number, log: boolean): Range =>
        log ? [c * Math.pow(r[0] / c, k), c * Math.pow(r[1] / c, k)] : [c + (r[0] - c) * k, c + (r[1] - c) * k];
      setView({ x: zoom(v.x, G.ix(p.x), mx.log), y: zoom(v.y, G.iy(p.y), my.log) });
    };
    el.addEventListener("wheel", h, { passive: false });
    return () => el.removeEventListener("wheel", h);
  });

  return (
    <div ref={wrap} className="fill scatter" style={{ minHeight: 340 }}>
      <canvas ref={cv} onMouseDown={onDown} onMouseMove={onMove} onMouseUp={onUp}
        onMouseLeave={() => { tip(null); if (drag.current) { drag.current = null; redraw((n) => n + 1); } }}
        onDoubleClick={() => setView(null)} />
      {view && <button className="zoomchip scatter-reset" onClick={() => setView(null)}>zoomed  ✕</button>}
      <div className="scatter-hint">drag to zoom · scroll to zoom · shift-drag to pan · double-click to reset · drag an orange line to move that threshold</div>
    </div>
  );
}

/* ---------- Spearman correlation matrix ---------- */
function Correlations({ met, vals, onPick }: { met: Metric[]; vals: Vals; onPick: (x: string, y: string) => void }) {
  const theme = useStore((s) => s.theme);
  const wrap = useRef<HTMLDivElement>(null), cv = useRef<HTMLCanvasElement>(null);
  const geo = useRef<{ left: number; top: number; cell: number; m: number } | null>(null);
  const W = useWidth(wrap);
  const rho = useMemo(() => met.map((a) => met.map((b) => (a === b ? 1 : spearman(vals[a.key], vals[b.key])))), [met, vals]);
  useEffect(() => {
    const m = met.length, left = 136, top = 8;
    const cell = Math.max(14, Math.min(34, Math.floor((W - left - 12) / Math.max(1, m)))), H = top + cell * m + 132;
    const g = fit(cv.current!, W, H);
    g.clearRect(0, 0, W, H);
    if (!m) return;
    const neg = rgbOf(DIVERGING[0]), pos = rgbOf("#E34A33"), mid = rgbOf("--panel2");
    geo.current = { left, top, cell, m };
    g.font = `11px ${FONT}`;
    for (let r = 0; r < m; r++)
      for (let c = 0; c < m; c++) {
        const v = rho[r][c], x = left + c * cell, y = top + r * cell;
        g.fillStyle = r === c ? css("--field") : Number.isNaN(v) ? css("--cnna") : v >= 0 ? mixRGB(mid, pos, v) : mixRGB(mid, neg, -v);
        g.beginPath();
        if (g.roundRect) g.roundRect(x + 1, y + 1, cell - 2, cell - 2, 3); else g.rect(x + 1, y + 1, cell - 2, cell - 2);
        g.fill();
        if (cell >= 28 && !Number.isNaN(v) && r !== c) {
          g.fillStyle = Math.abs(v) > 0.6 ? "#fff" : css("--ink2"); g.textAlign = "center"; g.textBaseline = "middle"; g.font = `10px ${FONT}`;
          g.fillText(v.toFixed(1).replace("-0.0", "0.0"), x + cell / 2, y + cell / 2);
          g.font = `11px ${FONT}`;
        }
      }
    g.fillStyle = css("--ink2"); g.textAlign = "right"; g.textBaseline = "middle";
    met.forEach((mt, r) => g.fillText(mt.label, left - 8, top + r * cell + cell / 2));
    met.forEach((mt, c) => { g.save(); g.translate(left + c * cell + cell / 2, top + m * cell + 6); g.rotate(-Math.PI / 4); g.fillText(mt.label, 0, 0); g.restore(); });
    const kx = left, ky = H - 12, kw = Math.min(220, cell * m);
    for (let k = 0; k < kw; k++) { const t = (k / (kw - 1)) * 2 - 1; g.fillStyle = t >= 0 ? mixRGB(mid, pos, t) : mixRGB(mid, neg, -t); g.fillRect(kx + k, ky, 1, 8); }
    g.fillStyle = css("--ink3"); g.textBaseline = "bottom";
    g.textAlign = "left"; g.fillText("−1", kx, ky - 2); g.textAlign = "center"; g.fillText("0", kx + kw / 2, ky - 2); g.textAlign = "right"; g.fillText("+1", kx + kw, ky - 2);
  }, [W, met, rho, theme]);
  const hit = (e: React.MouseEvent): [number, number] | null => {
    const c = geo.current;
    if (!c) return null;
    const r = cv.current!.getBoundingClientRect(), x = e.clientX - r.left - c.left, y = e.clientY - r.top - c.top;
    const i = Math.floor(y / c.cell), j = Math.floor(x / c.cell);
    return i >= 0 && j >= 0 && i < c.m && j < c.m && i !== j ? [i, j] : null;
  };
  const move = (e: React.MouseEvent) => {
    const h = hit(e);
    cv.current!.style.cursor = h ? "pointer" : "default";
    tip(e, h ? `<b>${esc(met[h[0]].label)}</b> × <b>${esc(met[h[1]].label)}</b><br><span class="k">Spearman ρ</span> ${rho[h[0]][h[1]].toFixed(2)}` : null);
  };
  return <div ref={wrap}><canvas ref={cv} onMouseMove={move} onMouseLeave={() => tip(null)} onClick={(e) => { const h = hit(e); if (h) onPick(met[h[1]].key, met[h[0]].key); }} /></div>;
}

/* ---------- one metric's histogram, passing and flagged stacked ---------- */
/** Arrow-key / spinner step of a filter's box: about 1/100 of the metric's range, a round number. */
function stepFor(key: string, range?: [number, number]): number {
  if (INT_METRICS.has(key)) return 1;
  const span = range ? range[1] - range[0] : 1;
  return Math.pow(10, Math.floor(Math.log10(span / 100 || 1)));
}

/** Where a new minimum or maximum starts: median ∓ 3 MAD, kept inside the data. */
function startBound(m: Metric, side: "min" | "max", range: [number, number]): number {
  const v = side === "max" ? m.hi : m.lo;
  const inside = Number.isFinite(v) && v > range[0] && v < range[1];
  const x = inside ? v : side === "max" ? range[0] + 0.9 * (range[1] - range[0]) : range[0] + 0.1 * (range[1] - range[0]);
  return INT_METRICS.has(m.key) ? Math.round(x) : +x.toPrecision(3);
}

/* ---------- one metric's histogram and its filter, shared by the card and the window ---------- */
interface HistGeo { X: (v: number) => number; inv: (x: number) => number }

/** Draw a metric's histogram: passing, flagged by another filter, flagged by this one. */
function drawHist(g: CanvasRenderingContext2D, W: number, H: number, m: Metric, values: Float64Array, passing: boolean[], f: Filter, big: boolean): HistGeo | null {
  g.clearRect(0, 0, W, H);
  const sc = Array.from(values, (v) => (m.log ? (v > 0 ? Math.log10(v) : NaN) : v)), ok = sc.filter((v) => !Number.isNaN(v));
  if (!ok.length) return null;
  // the axis spans the 0.5th–99.5th percentile: a few extreme cells pile into the end bars
  const sorted = [...ok].sort((x, y) => x - y), q = (p: number) => sorted[Math.round(p * (sorted.length - 1))];
  let lo = q(0.005), hi = q(0.995);
  if (!(hi > lo)) { lo = sorted[0]; hi = sorted[sorted.length - 1]; }
  if (lo === hi) { lo -= 0.5; hi += 0.5; }
  const clipLo = sorted[0] < lo, clipHi = sorted[sorted.length - 1] > hi;
  const fs = big ? 12.5 : 11, pad = big ? 14 : 8;
  const box = { l: pad, t: fs + 10, w: W - 2 * pad, h: H - (fs + 10) - (fs + 12) };
  const nb = big ? Math.min(80, Math.max(30, Math.round(box.w / 9))) : 26;
  const P = new Array(nb).fill(0), O = new Array(nb).fill(0), F = new Array(nb).fill(0);
  sc.forEach((v, i) => {
    if (Number.isNaN(v)) return;
    const k = Math.max(0, Math.min(nb - 1, Math.floor(((v - lo) / (hi - lo)) * nb)));
    (passing[i] ? P : fails(f, values[i]) ? F : O)[k]++;
  });
  const ymax = Math.max(...P.map((p, k) => p + O[k] + F[k])) || 1, bw = box.w / nb;
  const tr = (v: number) => (m.log ? Math.log10(v) : v), back = (v: number) => (m.log ? Math.pow(10, v) : v);
  const X = (v: number) => box.l + ((tr(v) - lo) / (hi - lo)) * box.w;
  // the excluded side(s)
  g.fillStyle = css("--qc-fail"); g.globalAlpha = 0.07;
  if (f.min != null) g.fillRect(box.l, box.t, Math.max(0, Math.min(box.w, X(f.min) - box.l)), box.h);
  if (f.max != null) { const x = Math.max(box.l, Math.min(box.l + box.w, X(f.max))); g.fillRect(x, box.t, box.l + box.w - x, box.h); }
  g.globalAlpha = 1;
  const gap = bw > 6 ? 2 : 1;
  for (let k = 0; k < nb; k++) {
    let y = box.t + box.h;
    const x = box.l + k * bw + gap / 2, w = Math.max(1, bw - gap);
    for (const [n, c] of [[P[k], "--qc-pass"], [O[k], "--qc-other"], [F[k], "--qc-fail"]] as const) {
      if (!n) continue;
      const h = box.h * (n / ymax);
      g.fillStyle = css(c);
      roundTop(g, x, y - h, w, h - (y < box.t + box.h ? 1 : 0), Math.min(2, w / 2));
      y -= h;
    }
  }
  // median: a thin line, labelled at the top
  const xm = X(m.median);
  g.strokeStyle = css("--ink"); g.globalAlpha = 0.45; g.lineWidth = 1; g.setLineDash([3, 3]);
  g.beginPath(); g.moveTo(xm, box.t); g.lineTo(xm, box.t + box.h); g.stroke();
  g.setLineDash([]); g.globalAlpha = 1;
  g.font = `${fs}px ${FONT}`; g.textBaseline = "alphabetic"; g.fillStyle = css("--ink3");
  const ml = `median ${fmtMetric(m.key, m.median)}`, mw = g.measureText(ml).width;
  g.textAlign = "left"; g.fillText(ml, Math.max(box.l, Math.min(box.l + box.w - mw, xm - mw / 2)), box.t - 5);
  // filter handles
  for (const v of [f.min, f.max]) {
    if (v == null) continue;
    const x = Math.max(box.l, Math.min(box.l + box.w, X(v)));
    g.strokeStyle = css("--qc-fail"); g.lineWidth = 2;
    g.beginPath(); g.moveTo(x, box.t + 4); g.lineTo(x, box.t + box.h); g.stroke();
    g.fillStyle = css("--qc-fail"); g.beginPath(); g.arc(x, box.t + 4, big ? 6 : 4.5, 0, 7); g.fill();
    g.fillStyle = css("--panel"); g.beginPath(); g.arc(x, box.t + 4, big ? 2.2 : 1.6, 0, 7); g.fill();
  }
  // axis ends, and the ticks between them in the big view
  const ly = H - 6;
  g.fillStyle = css("--ink3"); g.textBaseline = "alphabetic";
  const loL = (clipLo ? "< " : "") + fmtMetric(m.key, back(lo)), hiL = (clipHi ? "> " : "") + fmtMetric(m.key, back(hi));
  g.textAlign = "left"; g.fillText(loL, box.l, ly);
  g.textAlign = "right"; g.fillText(hiL, box.l + box.w, ly);
  if (big) {
    const loW = g.measureText(loL).width, hiW = g.measureText(hiL).width;
    g.textAlign = "center";
    for (let i = 1; i < 6; i++) {
      const t = lo + (i / 6) * (hi - lo), x = box.l + (i / 6) * box.w, lab = fmtMetric(m.key, back(t));
      const w = g.measureText(lab).width;
      if (x - w / 2 > box.l + loW + 10 && x + w / 2 < box.l + box.w - hiW - 10) g.fillText(lab, x, ly);
    }
  }
  return { X, inv: (x) => back(lo + ((x - box.l) / box.w) * (hi - lo)) };
}

/** Canvas, drawing and drag handling of a metric's histogram. */
function useFilterHist(m: Metric, values: Float64Array, passing: boolean[], f: Filter, H: number, big: boolean, onSide: (side: "min" | "max", v: number) => void) {
  const theme = useStore((s) => s.theme);
  const cv = useRef<HTMLCanvasElement>(null), wrap = useRef<HTMLDivElement>(null);
  const W = useWidth(wrap, 240);
  const geo = useRef<HistGeo | null>(null);
  const drag = useRef<"min" | "max" | null>(null);
  useEffect(() => {
    geo.current = drawHist(fit(cv.current!, W, H), W, H, m, values, passing, { min: f.min, max: f.max }, big);
  }, [W, H, m, values, passing, theme, f.min, f.max, big]);
  const xOf = (e: React.PointerEvent) => e.clientX - cv.current!.getBoundingClientRect().left;
  const near = (x: number): "min" | "max" | null => {
    const G = geo.current;
    if (!G) return null;
    for (const side of ["min", "max"] as const) { const v = f[side]; if (v != null && Math.abs(G.X(v) - x) <= (big ? 9 : 7)) return side; }
    return null;
  };
  const handlers = {
    onPointerDown: (e: React.PointerEvent) => {
      const s = near(xOf(e));
      if (!s) return;
      drag.current = s;
      try { (e.target as HTMLElement).setPointerCapture(e.pointerId); } catch { /* synthetic pointer */ }
    },
    onPointerMove: (e: React.PointerEvent) => {
      const x = xOf(e), G = geo.current;
      if (drag.current && G) onSide(drag.current, G.inv(x));
      else (e.target as HTMLElement).style.cursor = near(x) ? "ew-resize" : "default";
    },
    onPointerUp: () => { drag.current = null; },
  };
  return { cv, wrap, handlers };
}

function useFilterEdit(m: Metric, filter: Filter | undefined, range: [number, number] | undefined, onChange: (f: Filter | null) => void) {
  const f: Filter = filter ?? { min: null, max: null };
  const clampVal = (v: number) => {
    if (range) v = Math.max(range[0], Math.min(range[1], v));
    return INT_METRICS.has(m.key) ? Math.round(v) : +v.toPrecision(3);
  };
  const setSide = (side: "min" | "max", v: number | null) => {
    const nf = { ...f, [side]: v == null ? null : clampVal(v) };
    if (nf.min != null && nf.max != null && nf.min > nf.max) return;
    onChange(nf.min == null && nf.max == null ? null : nf);
  };
  const bound = (side: "min" | "max") => {
    const v = f[side];
    return (
      <div className={"fbound" + (v != null ? " on" : "")}>
        <button className="fbtn" disabled={!range} title={v != null ? `Remove the ${side}imum` : `Add a ${side}imum`}
          onClick={() => setSide(side, v != null ? null : startBound(m, side, range!))}>{side === "min" ? "Min" : "Max"}</button>
        {v != null && (
          <input className="fnum" type="number" value={v} step={stepFor(m.key, range)} aria-label={`${m.label} ${side}imum`}
            onChange={(e) => { const x = parseFloat(e.target.value); setSide(side, Number.isFinite(x) ? x : null); }} />
        )}
      </div>
    );
  };
  return { f, setSide, bound };
}

interface FilterProps {
  m: Metric; values: Float64Array; passing: boolean[]; filter?: Filter; count: number; range?: [number, number];
  onChange: (f: Filter | null) => void;
}

/* ---------- a metric's card in the grid; double-click opens it large ---------- */
export function FilterCard({ plotted, onPlot, onOpen, ...p }: FilterProps & { plotted: boolean; onPlot: () => void; onOpen: () => void }) {
  const { f, setSide, bound } = useFilterEdit(p.m, p.filter, p.range, p.onChange);
  const { cv, wrap, handlers } = useFilterHist(p.m, p.values, p.passing, f, 128, false, setSide);
  return (
    <div id={`filter-${p.m.key}`} className={"dist" + (isActive(f) ? " filtered" : "") + (plotted ? " on" : "")} onDoubleClick={(e) => { if (!(e.target as HTMLElement).closest("input, button")) onOpen(); }} title="Double-click to open it large">
      <div className="dhead">
        <button className="dn" title={p.m.help + " Click to plot it in Explore."} onClick={onPlot}>{p.m.label}</button>
        {isActive(f) ? <span className={"fcount" + (p.count ? " bad" : "")}>{p.count.toLocaleString()} flagged</span> : <span className="dg">{p.m.group}</span>}
      </div>
      <div ref={wrap}><canvas ref={cv} {...handlers} /></div>
      <div className="dfoot">
        {bound("min")}{bound("max")}
        <button className="fopen" onClick={onOpen} aria-label={`Open ${p.m.label} large`} title="Open large">⤢</button>
      </div>
    </div>
  );
}

/* ---------- the large window of one metric ---------- */
export function FilterWindow({ onClose, onPrev, onNext, onPlot, ...p }: FilterProps & { onClose: () => void; onPrev: () => void; onNext: () => void; onPlot: () => void }) {
  const { f, setSide, bound } = useFilterEdit(p.m, p.filter, p.range, p.onChange);
  const { cv, wrap, handlers } = useFilterHist(p.m, p.values, p.passing, f, 360, true, setSide);
  useEffect(() => {
    const h = (e: KeyboardEvent) => {
      if ((e.target as HTMLElement).tagName === "INPUT") return;
      if (e.key === "Escape") { e.stopPropagation(); onClose(); }
      else if (e.key === "ArrowLeft") onPrev();
      else if (e.key === "ArrowRight") onNext();
    };
    addEventListener("keydown", h, true);
    return () => removeEventListener("keydown", h, true);
  }, [onClose, onPrev, onNext]);
  const m = p.m, stat = (k: string, v: string) => <div className="wstat"><span>{k}</span><b>{v}</b></div>;
  return (
    <div className="modal-back" onMouseDown={(e) => e.target === e.currentTarget && onClose()}>
      <div className="fwin" role="dialog" aria-label={`${m.label} filter`}>
        <div className="fwin-head">
          <button className="iconbtn" onClick={onPrev} aria-label="Previous metric">‹</button>
          <div className="fwin-title"><h2>{m.label}</h2><span className="ok">{m.group} · {m.dir === "upper" ? "bad when high" : m.dir === "lower" ? "bad when low" : m.dir === "none" ? "describes the cell, not its quality" : "bad on either side"}</span></div>
          <button className="iconbtn" onClick={onNext} aria-label="Next metric">›</button>
          <span className="spacer" />
          {isActive(f) && <span className={"fcount big" + (p.count ? " bad" : "")}>{p.count.toLocaleString()} cells flagged</span>}
          <button className="iconbtn" onClick={onPlot}>Plot in Explore</button>
          <button className="iconbtn" onClick={onClose} aria-label="Close">✕</button>
        </div>
        <p className="fwin-help">{m.help}</p>
        <div ref={wrap} className="fwin-plot"><canvas ref={cv} {...handlers} /></div>
        <div className="fwin-foot">
          <div className="fwin-ctl">{bound("min")}{bound("max")}</div>
          <div className="qc-key"><span><i style={{ background: "var(--qc-pass)" }} />passes</span><span><i style={{ background: "var(--qc-fail)" }} />flagged by this filter</span><span><i style={{ background: "var(--qc-other)" }} />flagged by another</span></div>
          <div className="wstats">
            {stat("cells", m.n.toLocaleString())}{stat("median", fmtMetric(m.key, m.median))}{stat("MAD", m.mad == null ? "log" : fmtMetric(m.key, m.mad))}
            {stat("min", fmtMetric(m.key, m.min))}{stat("max", fmtMetric(m.key, m.max))}
          </div>
        </div>
        <p className="hint">Drag the handle, or type a value · ← → other metrics · Esc closes</p>
      </div>
    </div>
  );
}

/** The Cell QC tab: flags, an explorer of any two metrics, their correlations, the distribution of each and their statistics. */
export function QCTab({ children }: { children: ReactNode }) {
  const d = useDerived();
  const { meta, vals, metAll, passing, reasons, shown, ranges, failCounts, nPass, filters } = d;
  const selected = useStore((s) => s.selected);
  const set = useStore((s) => s.set);
  const setFilter = useStore((s) => s.setFilter);
  const shownSet = useMemo(() => new Set(shown), [shown]);
  const met = useMemo(() => metAll.filter((m) => new Set(Array.from(vals[m.key]).filter((v) => !Number.isNaN(v))).size > 1), [metAll, vals]);
  const [ax, setAx] = useState<{ x: string | null; y: string | null }>({ x: null, y: null });
  const [openKey, setOpenKey] = useState<string | null>(null);
  const has = (k: string) => met.some((m) => m.key === k);
  const xKey = ax.x && has(ax.x) ? ax.x : has("reads_per_bin") ? "reads_per_bin" : has("mapd") ? "mapd" : met[0]?.key;
  const yKey = ax.y && has(ax.y) ? ax.y : has("mapd") && xKey !== "mapd" ? "mapd" : has("median_residual") ? "median_residual" : met[1]?.key ?? xKey;
  const mx = met.find((m) => m.key === xKey), my = met.find((m) => m.key === yKey);
  const groups = [...new Set(metAll.map((m) => m.group))];
  const dirTxt = { upper: "bad when high", lower: "bad when low", both: "bad on either side", none: "describes the cell" };

  const tiles: [string, string, string][] = [["Cells passing", nPass.toLocaleString(), `of ${meta.n_cells.toLocaleString()} (${meta.n_cells ? Math.round((100 * nPass) / meta.n_cells) : 0}%)`]];
  for (const [k, f] of Object.entries(filters)) if (isActive(f) && ranges[k]) tiles.push([metAll.find((m) => m.key === k)?.label ?? k, (failCounts[k] ?? 0).toLocaleString(), `flagged · ${describe(k, f)}`]);
  const nActive = Object.values(filters).filter(isActive).length;

  if (!metAll.length) return <><div className="card"><div className="empty">This dataset has no per-cell QC metrics.</div></div>{children}</>;
  return (
    <>
      <div className="tiles">
        {tiles.map(([k, v, s]) => <div className="tile" key={k}><div className="k">{k}</div><div className="v">{v}</div><div className="s">{s}</div></div>)}
      </div>
      {(() => {
        const i = met.findIndex((m) => m.key === openKey);
        if (i < 0) return null;
        const m = met[i], go = (d: number) => setOpenKey(met[(i + d + met.length) % met.length].key);
        return (
          <FilterWindow m={m} values={vals[m.key]} passing={passing} filter={filters[m.key]} count={failCounts[m.key] ?? 0} range={ranges[m.key]}
            onChange={(f) => setFilter(m.key, f)} onClose={() => setOpenKey(null)} onPrev={() => go(-1)} onNext={() => go(1)}
            onPlot={() => { setAx((a) => ({ ...a, x: m.key })); setOpenKey(null); }} />
        );
      })()}
      <div className="card" id="qc-filters">
        <div className="hd">
          <h2>Filters</h2>
          <span className="sub">{nPass.toLocaleString()} of {meta.n_cells.toLocaleString()} cells pass · drag a handle, switch Min / Max on, double-click a histogram to open it large</span>
          <span className="qc-key"><span><i style={{ background: "var(--qc-pass)" }} />passes</span><span><i style={{ background: "var(--qc-fail)" }} />flagged by this filter</span><span><i style={{ background: "var(--qc-other)" }} />flagged by another</span></span>
          <div className="fbar">
            <button className="iconbtn" onClick={() => set({ filters: defaultFilters(vals) })} title="The ASCAT.sc Explorer's defaults: median + 2 MAD for the residual and MAPD, percentiles for coverage, CNA bins and Gini">Defaults</button>
            <button className="iconbtn" disabled={!nActive} onClick={() => set({ filters: {} })}>Clear all</button>
            <button className="iconbtn" onClick={() => downloadText(cellsTable(meta, vals, filters, d), `${(meta.name || "cnvscope").replace(/[^\w.-]+/g, "_")}_cells_qc.tsv`)} title="Every cell with pass/fail, the failed filters and all metrics (TSV)"><DlIcon />Download cells</button>
          </div>
        </div>
        <div className="bd">
        <div className="segrow"><SegmentControl compact /></div>
        <div className="dists">
          {met.map((m) => (
            <FilterCard key={m.key} m={m} values={vals[m.key]} passing={passing} filter={filters[m.key]} count={failCounts[m.key] ?? 0}
              range={ranges[m.key]} plotted={m.key === xKey} onPlot={() => setAx((a) => ({ ...a, x: m.key }))}
              onOpen={() => setOpenKey(m.key)} onChange={(f) => setFilter(m.key, f)} />
          ))}
        </div>
        </div>
      </div>
      <div className="grid2">
        <div className="card">
          <div className="hd">
            <h2>Explore</h2><span className="sub">each dot is a cell · click to open it</span>
            <div className="axsel">
              <label>X <select aria-label="x axis metric" value={xKey || ""} onChange={(e) => setAx((a) => ({ ...a, x: e.target.value }))}>{met.map((m) => <option key={m.key} value={m.key}>{m.label}</option>)}</select></label>
              <label>Y <select aria-label="y axis metric" value={yKey || ""} onChange={(e) => setAx((a) => ({ ...a, y: e.target.value }))}>{met.map((m) => <option key={m.key} value={m.key}>{m.label}</option>)}</select></label>
            </div>
          </div>
          <div className="bd">
            <Scatter cells={meta.cells} mx={mx} my={my} vals={vals} passing={passing} shown={shownSet} reasons={reasons} sel={selected} onSelect={(i) => set({ selected: i })} />
          </div>
        </div>
        <div className="card">
          <div className="hd"><h2>Metric correlations</h2><span className="sub">Spearman ρ · click a square to plot the pair</span></div>
          <div className="bd"><Correlations met={met} vals={vals} onPick={(x, y) => setAx({ x, y })} /></div>
        </div>
      </div>
      {children}
      <div className="card">
        <div className="hd"><h2>Statistics</h2><span className="sub">per metric, all cells</span></div>
        <div className="bd tablewrap"><table>
          <thead><tr><th>metric</th><th>cells</th><th>median</th><th>MAD</th><th>min</th><th>max</th><th title="cells beyond median ± 3 MAD on the bad side (log scale for reads)">outliers</th></tr></thead>
          <tbody>
            {groups.flatMap((gname) => [
              <tr className="grp" key={gname}><td colSpan={7}>{gname}</td></tr>,
              ...metAll.filter((m) => m.group === gname).map((m) => (
                <tr key={m.key}>
                  <td title={m.help}>{m.label} <span className="ok">· {dirTxt[m.dir]}</span></td>
                  <td>{m.n.toLocaleString()}</td><td>{fmtMetric(m.key, m.median)}</td>
                  <td>{m.log ? <span className="ok">log</span> : fmtMetric(m.key, m.mad)}</td>
                  <td>{fmtMetric(m.key, m.min)}</td><td>{fmtMetric(m.key, m.max)}</td>
                  <td>{m.out ? <><span style={{ color: "var(--serious)", fontWeight: 600 }}>{m.out}</span> <span className="ok">({((100 * m.out) / m.n).toFixed(1)}%)</span></> : <span className="ok">0</span>}</td>
                </tr>
              )),
            ])}
          </tbody>
        </table></div>
      </div>
    </>
  );
}
