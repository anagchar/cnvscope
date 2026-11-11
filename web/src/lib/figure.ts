// Publication figures of the heatmap and of a cell's profile, drawn on a 2D canvas in
// points (1/72 inch) so the same code renders the dialog's preview and the 300/600 dpi file.

import type { CellDetail, Meta } from "../api/client";
import { chrLabel, css, niceTicks } from "./canvas";
import { layoutVisible } from "./dendrogram";
import { STAIN_COLORS, binAt, binColumns, binRange, chromLabels, type Band } from "./genome";
import { buildLut, discreteColor, valueToUnit, type Scale } from "./scales";

export const FIG_FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif";

export interface FigTheme {
  bg: string; ink: string; ink2: string; ink3: string; grid: string; sep: string; frame: string;
  na: string; mean: string; call: string; major: string; minor: string; points: string; bars: string;
}

/** "print": white background and dark text for journals; "screen": the app's current theme. */
export function figTheme(print: boolean): FigTheme {
  if (print)
    return {
      bg: "#FFFFFF", ink: "#1A1A1E", ink2: "#4A4A52", ink3: "#808089", grid: "#E8E8EC", sep: "#FFFFFF", frame: "#B9B9C1",
      na: "#F3F3F5", mean: "#2B2B30", call: "#2F2F35", major: "#D55E00", minor: "#0072B2", points: "#8A8A92", bars: "#5E7FA8",
    };
  return {
    bg: css("--panel"), ink: css("--ink"), ink2: css("--ink2"), ink3: css("--ink3"), grid: css("--grid"), sep: css("--hm-gap"), frame: css("--line2"),
    na: css("--cnna"), mean: css("--mean"), call: css("--call"), major: css("--major"), minor: css("--minor"), points: css("--points"), bars: css("--series1"),
  };
}

export interface CommonOpts {
  /** figure size in points */
  w: number;
  h: number;
  fontPt: number;
  title: string;
  subtitle: string;
  theme: FigTheme;
  legend: boolean;
}

const rgb = (s: string): [number, number, number] => {
  const c = document.createElement("canvas").getContext("2d")!;
  c.fillStyle = s;
  c.fillRect(0, 0, 1, 1);
  const d = c.getImageData(0, 0, 1, 1).data;
  return [d[0], d[1], d[2]];
};

function header(g: CanvasRenderingContext2D, o: CommonOpts, pad: number): number {
  let y = pad;
  if (o.title) {
    g.fillStyle = o.theme.ink; g.font = `600 ${o.fontPt * 1.25}px ${FIG_FONT}`; g.textAlign = "left"; g.textBaseline = "top";
    g.fillText(o.title, pad, y);
    y += o.fontPt * 1.25 * 1.3;
  }
  if (o.subtitle) {
    g.fillStyle = o.theme.ink3; g.font = `${o.fontPt * 0.9}px ${FIG_FONT}`; g.textBaseline = "top";
    g.fillText(o.subtitle, pad, y);
    y += o.fontPt * 1.4;
  }
  return y === pad ? pad : y + o.fontPt * 0.4;
}

/** A row of colour swatches with labels under them; returns the width used. */
function swatches(g: CanvasRenderingContext2D, x: number, y: number, fs: number, items: [string, string][], theme: FigTheme, maxW: number): number {
  const sw = fs * 1.7, sh = fs * 0.85, gap = fs * 0.35;
  g.font = `${fs * 0.85}px ${FIG_FONT}`; g.textAlign = "center"; g.textBaseline = "top";
  let cx = x;
  for (const [col, label] of items) {
    const w = Math.max(sw, g.measureText(label).width);
    if (cx + w > x + maxW) break;
    g.fillStyle = col; g.fillRect(cx + (w - sw) / 2, y, sw, sh);
    g.strokeStyle = theme.frame; g.lineWidth = 0.3; g.strokeRect(cx + (w - sw) / 2, y, sw, sh);
    g.fillStyle = theme.ink2; g.fillText(label, cx + w / 2, y + sh + fs * 0.2);
    cx += w + gap;
  }
  return cx - x;
}

function scaleLegend(g: CanvasRenderingContext2D, x: number, y: number, fs: number, title: string, s: Scale, theme: FigTheme, maxW: number) {
  g.fillStyle = theme.ink; g.font = `600 ${fs * 0.9}px ${FIG_FONT}`; g.textAlign = "left"; g.textBaseline = "top";
  g.fillText(title, x, y);
  const tw = g.measureText(title).width + fs;
  if (s.kind === "discrete") {
    swatches(g, x + tw, y - fs * 0.05, fs, [...s.stops.map((c, i) => [c.startsWith("--") ? css(c) : c, s.labels[i]] as [string, string]), [theme.na, "NA"]], theme, maxW - tw);
  } else {
    const bw = Math.min(fs * 12, maxW - tw), bh = fs * 0.85;
    const grad = g.createLinearGradient(x + tw, 0, x + tw + bw, 0);
    s.stops.forEach((c, i) => grad.addColorStop(i / (s.stops.length - 1), c.startsWith("--") ? css(c) : c));
    g.fillStyle = grad; g.fillRect(x + tw, y, bw, bh);
    g.fillStyle = theme.ink2; g.font = `${fs * 0.85}px ${FIG_FONT}`; g.textBaseline = "top";
    g.textAlign = "left"; g.fillText(s.labels[0], x + tw, y + bh + fs * 0.2);
    g.textAlign = "right"; g.fillText(s.labels[s.labels.length - 1], x + tw + bw, y + bh + fs * 0.2);
  }
}

/* ───────────────────────────── heatmap ───────────────────────────── */

export interface HeatmapFigure extends CommonOpts {
  meta: Meta;
  matrix: Float32Array;
  scale: Scale;
  scaleTitle: string;
  /** dataset rows shown, in display order */
  shown: number[];
  rowOf: Int32Array;
  /** genomic window and display rows to draw */
  gx: [number, number];
  rows: [number, number];
  dendro: { merge: [number, number][]; height: number[] } | null;
  cytobands: Band[] | null;
  annotation: { name: string; values: (string | number | boolean | null)[]; colors: Map<string, string> } | null;
  /** device pixels per point, for the matrix image */
  pxPerPt: number;
}

export function drawHeatmapFigure(g: CanvasRenderingContext2D, o: HeatmapFigure) {
  const { theme: T, fontPt: fs, meta } = o;
  const B = meta.bins, nb = meta.n_bins;
  const pad = fs * 1.2;
  g.fillStyle = T.bg; g.fillRect(0, 0, o.w, o.h);
  let top = header(g, o, pad);

  const [gx0, gx1] = o.gx, [r0, r1] = o.rows;
  const legendH = o.legend ? fs * 2.6 * (o.annotation ? 2 : 1) + fs * 0.6 : 0;
  const labelH = fs * 1.3, ideoH = o.cytobands ? fs * 0.75 : 0;
  const dendroW = o.dendro ? Math.min(o.w * 0.12, fs * 7) : 0;
  const annotW = o.annotation ? fs * 0.9 : 0;
  const x0 = pad + dendroW + (o.dendro ? fs * 0.2 : 0) + (annotW ? annotW + fs * 0.3 : 0);
  const hw = o.w - x0 - pad;
  top += labelH + (ideoH ? ideoH + fs * 0.35 : 0);
  const hh = o.h - top - legendH - pad;
  if (hw <= 10 || hh <= 10) return;
  const xs = (gx: number) => x0 + ((gx - gx0) / (gx1 - gx0)) * hw;
  const ys = (r: number) => top + ((r - r0) / (r1 - r0)) * hh;

  // matrix: one image pixel per device pixel, nearest bin and row
  const wpx = Math.max(1, Math.round(hw * o.pxPerPt)), hpx = Math.max(1, Math.round(hh * o.pxPerPt));
  const img = new ImageData(wpx, hpx);
  const lut = buildLut(o.scale, (t) => rgb(t.startsWith("--") ? css(t) : t));
  const na = rgb(T.na), sepc = rgb(T.sep);
  const disc = o.scale.kind === "discrete", n = disc ? lut.length / 4 : 256;
  const cols = binColumns(B, gx0, gx1, wpx, new Int32Array(wpx));
  const d = img.data;
  for (let py = 0; py < hpx; py++) {
    const row = o.shown[Math.min(o.shown.length - 1, Math.floor(r0 + ((py + 0.5) / hpx) * (r1 - r0)))];
    const base = row * nb;
    for (let px = 0; px < wpx; px++) {
      const b = cols[px], k = (py * wpx + px) * 4;
      let c: ArrayLike<number> = sepc, off = 0;
      if (b >= 0) {
        const v = o.matrix[base + b];
        if (Number.isNaN(v)) c = na;
        else {
          const u = valueToUnit(o.scale, v);
          c = lut; off = (disc ? u : Math.round(u * (n - 1))) * 4;
        }
      }
      d[k] = c[off]; d[k + 1] = c[off + 1]; d[k + 2] = c[off + 2]; d[k + 3] = 255;
    }
  }
  const off = document.createElement("canvas");
  off.width = wpx; off.height = hpx;
  off.getContext("2d")!.putImageData(img, 0, 0);
  g.imageSmoothingEnabled = false;
  g.drawImage(off, x0, top, hw, hh);
  g.imageSmoothingEnabled = true;

  // chromosome separators, labels, ideogram
  g.fillStyle = T.sep;
  for (const c of meta.chr_info) if (c.start_cum > gx0 && c.start_cum < gx1) g.fillRect(xs(c.start_cum) - 0.35, top, 0.7, hh);
  g.font = `${fs * 0.85}px ${FIG_FONT}`; g.fillStyle = T.ink2; g.textAlign = "center"; g.textBaseline = "alphabetic";
  const spans = meta.chr_info.filter((c) => c.end_cum >= gx0 && c.start_cum <= gx1).map((c) => [c.chr, xs(Math.max(c.start_cum, gx0)), xs(Math.min(c.end_cum, gx1))] as [string, number, number]);
  const labelY = top - (ideoH ? ideoH + fs * 0.35 : 0) - fs * 0.35;
  for (const [t, x] of chromLabels(g, spans)) g.fillText(t, x, labelY);
  if (o.cytobands && ideoH) {
    const iy = top - ideoH - fs * 0.2;
    const offs: Record<string, number> = {};
    meta.chr_info.forEach((c) => (offs[c.chr] = c.start_cum));
    for (const band of o.cytobands) {
      const ofs = offs[band.chr];
      if (ofs === undefined) continue;
      const s = Math.max(ofs + band.start, gx0), e = Math.min(ofs + band.end, gx1);
      if (e <= s) continue;
      g.fillStyle = STAIN_COLORS[band.stain] || "#AAAAAA";
      g.fillRect(xs(s), iy, Math.max(0.2, xs(e) - xs(s)), ideoH);
    }
    g.strokeStyle = T.frame; g.lineWidth = 0.3;
    for (const [, a, b] of spans) g.strokeRect(a, iy, b - a, ideoH);
  }
  g.strokeStyle = T.frame; g.lineWidth = 0.4; g.strokeRect(x0, top, hw, hh);

  // annotation strip
  if (o.annotation) {
    const ax = x0 - annotW - fs * 0.3;
    const vals = o.annotation.values;
    for (let r = Math.max(0, Math.floor(r0)); r < Math.min(o.shown.length, Math.ceil(r1)); r++) {
      const v = vals[o.shown[r]];
      if (v == null) continue;
      g.fillStyle = o.annotation.colors.get(String(v)) ?? T.na;
      const a = Math.max(top, ys(r)), b = Math.min(top + hh, ys(r + 1));
      g.fillRect(ax, a, annotW, Math.max(0.3, b - a + 0.15));
    }
  }

  // dendrogram
  if (o.dendro) {
    const lay = layoutVisible(o.dendro.merge, o.dendro.height, o.rowOf, r0, r1, o.shown.length);
    if (lay.maxH) {
      const dx = (h: number) => pad + dendroW - Math.sqrt(h / lay.maxH) * (dendroW - 1);
      g.strokeStyle = T.ink3; g.lineWidth = 0.45; g.beginPath();
      for (const s of lay.segments) {
        if (Math.max(s.y1, s.y2) < r0 || Math.min(s.y1, s.y2) > r1) continue;
        const a = s.y1 === s.y2 ? s.y1 : Math.max(r0, Math.min(r1, s.y1)), b = s.y1 === s.y2 ? s.y2 : Math.max(r0, Math.min(r1, s.y2));
        g.moveTo(dx(s.x1), ys(a)); g.lineTo(dx(s.x2), ys(b));
      }
      g.stroke();
    }
  }

  // legend
  if (o.legend) {
    let ly = top + hh + fs * 0.9;
    scaleLegend(g, x0, ly, fs, o.scaleTitle, o.scale, T, hw);
    if (o.annotation) {
      ly += fs * 2.6;
      g.fillStyle = T.ink; g.font = `600 ${fs * 0.9}px ${FIG_FONT}`; g.textAlign = "left"; g.textBaseline = "top";
      g.fillText(o.annotation.name, x0, ly);
      const tw = g.measureText(o.annotation.name).width + fs;
      swatches(g, x0 + tw, ly - fs * 0.05, fs, [...o.annotation.colors].map(([k, c]) => [c, k] as [string, string]), T, hw - tw);
    }
  }
}

/* ───────────────────────────── profile ───────────────────────────── */

export interface ProfileFigure extends CommonOpts {
  meta: Meta;
  cell: CellDetail;
  gx: [number, number];
  /** the integer call layer and the per-bin layer on the copy-number scale (already transformed) */
  calls: (number | null)[] | null;
  bins: (number | null)[] | null;
  allele: { major: (number | null)[]; minor: (number | null)[]; lo?: [(number | null)[], (number | null)[]]; hi?: [(number | null)[], (number | null)[]] } | null;
  reads: (number | null)[] | null;
  stateScale: Scale;
  /** colour bins by their state */
  colorBins: boolean;
}

type V = (number | null)[];
function runs(chr: string[], v: V): [number, number, number][] {
  const out: [number, number, number][] = [];
  for (let s = 0; s < v.length; ) {
    let e = s + 1;
    while (e < v.length && v[e] === v[s] && chr[e] === chr[s]) e++;
    if (v[s] != null) out.push([s, e, v[s] as number]);
    s = e;
  }
  return out;
}

export function drawProfileFigure(g: CanvasRenderingContext2D, o: ProfileFigure) {
  const { theme: T, fontPt: fs, meta } = o;
  const B = meta.bins;
  const pad = fs * 1.2;
  g.fillStyle = T.bg; g.fillRect(0, 0, o.w, o.h);
  const top = header(g, o, pad);
  const legendH = o.legend ? fs * 2 : 0;
  const readsH = o.reads ? (o.h - top) * 0.24 : 0;
  const left = pad + fs * 2.6, right = o.w - pad;
  const plotH = o.h - top - legendH - readsH - pad - fs * 1.6 - (readsH ? fs * 0.8 : 0);
  if (plotH < 10) return;
  const [gx0, gx1] = o.gx;
  const [b0, b1] = binRange(B.start_cum, B.end_cum, gx0, gx1);
  const xs = (gx: number) => left + ((gx - gx0) / (gx1 - gx0)) * (right - left);
  const mx = (v: V | null | undefined) => { let m = 0; if (v) for (let i = b0; i < b1; i++) { const x = v[i]; if (x != null && x > m) m = x; } return m; };
  const ymax = Math.max(5, Math.min(12, (o.allele ? mx(o.allele.major) : mx(o.calls)) + 1));
  const Y = (v: number) => top + plotH - (Math.max(-0.3, Math.min(ymax + 0.4, v)) / (ymax + 0.4)) * plotH;

  // grid, axis
  g.font = `${fs * 0.85}px ${FIG_FONT}`; g.textAlign = "right"; g.textBaseline = "middle"; g.lineWidth = 0.4;
  for (let y = 0; y <= ymax; y++) {
    g.strokeStyle = y === 2 ? T.ink3 : T.grid; g.setLineDash(y === 2 ? [] : [1.5, 1.5]);
    g.beginPath(); g.moveTo(left, Y(y)); g.lineTo(right, Y(y)); g.stroke();
    g.fillStyle = T.ink2; g.fillText(String(y), left - fs * 0.4, Y(y));
  }
  g.setLineDash([]);
  g.save(); g.translate(pad + fs * 0.5, top + plotH / 2); g.rotate(-Math.PI / 2); g.textAlign = "center"; g.fillStyle = T.ink2; g.fillText("Copy number", 0, 0); g.restore();
  const spans: [string, number, number][] = [];
  g.strokeStyle = T.frame; g.lineWidth = 0.4;
  for (const c of meta.chr_info) {
    if (c.end_cum < gx0 || c.start_cum > gx1) continue;
    const a = xs(Math.max(c.start_cum, gx0)), b = xs(Math.min(c.end_cum, gx1));
    if (c.start_cum > gx0) { g.beginPath(); g.moveTo(a, top); g.lineTo(a, top + plotH); g.stroke(); }
    spans.push([c.chr, a, b]);
  }
  g.fillStyle = T.ink2; g.textAlign = "center"; g.textBaseline = "top";
  for (const [t, x] of chromLabels(g, spans)) g.fillText(t, x, top + plotH + fs * 0.35);

  g.save(); g.beginPath(); g.rect(left, top, right - left, plotH); g.clip();
  const nz = b1 - b0;
  if (o.bins) {
    const r = Math.max(0.35, Math.min(1.4, ((right - left) / Math.max(1, nz)) * 0.45));
    g.globalAlpha = 0.85;
    for (let i = b0; i < b1; i++) {
      const v = o.bins[i];
      if (v == null || !Number.isFinite(v)) continue;
      g.fillStyle = (o.colorBins && o.calls && discreteColor(o.stateScale, o.calls[i] ?? NaN)) || T.points;
      g.beginPath(); g.arc(xs((B.start_cum[i] + B.end_cum[i]) / 2), Y(v), r, 0, 7); g.fill();
    }
    g.globalAlpha = 1;
  }
  const seg = (v: V, col: string, w: number, off = 0) => {
    g.strokeStyle = col; g.lineWidth = w; g.lineCap = "round"; g.beginPath();
    for (const [s, e, x] of runs(B.chr, v)) { if (e <= b0 || s >= b1) continue; g.moveTo(xs(B.start_cum[s]), Y(Math.min(x, ymax) + off)); g.lineTo(xs(B.end_cum[e - 1]), Y(Math.min(x, ymax) + off)); }
    g.stroke();
  };
  if (o.allele) {
    const A = o.allele;
    for (const [arr, col, k] of [[A.major, T.major, 0], [A.minor, T.minor, 1]] as const) {
      if (A.lo && A.hi) {
        g.fillStyle = col; g.globalAlpha = 0.18;
        for (let i = b0; i < b1; i++) {
          const a = arr[i], l = A.lo[k][i], h = A.hi[k][i];
          if (a == null || l == null || h == null) continue;
          const yt = Y(Math.max(a, h)), yb = Y(Math.min(a, l));
          g.fillRect(xs(B.start_cum[i]), yt, xs(B.end_cum[i]) - xs(B.start_cum[i]), yb - yt);
        }
        g.globalAlpha = 1;
      }
    }
    seg(A.major, T.major, fs * 0.22, 0.07);
    seg(A.minor, T.minor, fs * 0.22, -0.07);
  } else if (o.calls) seg(o.calls, T.call, fs * 0.2);
  g.restore();

  // read depth
  if (o.reads && readsH) {
    const rt = top + plotH + fs * 1.6 + fs * 0.8, rh = readsH - fs * 0.4;
    const seen: number[] = [];
    for (let i = b0; i < b1; i++) { const v = o.reads[i]; if (v != null) seen.push(v); }
    seen.sort((a, b) => a - b);
    const max = Math.max(1, seen.length ? seen[Math.floor(0.995 * (seen.length - 1))] * 1.1 : 1);
    const RY = (v: number) => rt + rh - (Math.min(v, max) / max) * rh;
    g.font = `${fs * 0.8}px ${FIG_FONT}`; g.textAlign = "right"; g.textBaseline = "middle"; g.lineWidth = 0.4;
    for (const t of niceTicks(0, max, 3)) { g.strokeStyle = T.grid; g.beginPath(); g.moveTo(left, RY(t)); g.lineTo(right, RY(t)); g.stroke(); g.fillStyle = T.ink3; g.fillText(String(t), left - fs * 0.4, RY(t)); }
    g.save(); g.translate(pad + fs * 0.5, rt + rh / 2); g.rotate(-Math.PI / 2); g.textAlign = "center"; g.fillStyle = T.ink2; g.fillText("Reads", 0, 0); g.restore();
    g.save(); g.beginPath(); g.rect(left, rt, right - left, rh); g.clip();
    g.fillStyle = T.bars; g.globalAlpha = 0.85;
    for (let i = b0; i < b1; i++) {
      const v = o.reads[i];
      if (v == null) continue;
      const a = xs(B.start_cum[i]), b = xs(B.end_cum[i]), y = RY(v);
      g.fillRect(a, y, Math.max(0.2, b - a), rt + rh - y);
    }
    g.restore();
  }

  // legend
  if (o.legend) {
    const ly = o.h - pad - fs * 1.1;
    g.font = `${fs * 0.85}px ${FIG_FONT}`; g.textAlign = "left"; g.textBaseline = "middle";
    let x = left;
    const item = (draw: () => void, label: string) => { draw(); g.fillStyle = T.ink2; g.fillText(label, x + fs * 1.6, ly); x += fs * 1.6 + g.measureText(label).width + fs * 1.2; };
    if (o.bins) item(() => {
      const cols = o.colorBins ? [o.stateScale.stops[0], o.stateScale.stops[2], o.stateScale.stops[4], o.stateScale.stops[6]] : [T.points];
      cols.forEach((c, i) => { g.fillStyle = c.startsWith("--") ? css(c) : c; g.beginPath(); g.arc(x + fs * 0.25 + i * fs * 0.35, ly, fs * 0.17, 0, 7); g.fill(); });
    }, o.colorBins ? "bins, coloured by state" : "bins");
    if (o.allele) {
      item(() => { g.fillStyle = T.major; g.fillRect(x, ly - fs * 0.11, fs * 1.2, fs * 0.22); }, "major allele");
      item(() => { g.fillStyle = T.minor; g.fillRect(x, ly - fs * 0.11, fs * 1.2, fs * 0.22); }, "minor allele");
      if (o.allele.lo) item(() => { g.fillStyle = T.major; g.globalAlpha = 0.25; g.fillRect(x, ly - fs * 0.3, fs * 1.2, fs * 0.6); g.globalAlpha = 1; }, "BAF 5–95% interval");
    } else if (o.calls) item(() => { g.fillStyle = T.call; g.fillRect(x, ly - fs * 0.1, fs * 1.2, fs * 0.2); }, "copy-number call");
  }
}

/** Genomic window label, e.g. "chr3:168.3 Mb – chr5:9.1 Mb". */
export function regionLabel(meta: Meta, gx: [number, number]): string {
  const B = meta.bins, a = binAt(B.end_cum, gx[0]), b = binAt(B.end_cum, gx[1]);
  const mb = (v: number) => (v / 1e6).toFixed(1) + " Mb";
  return `${chrLabel(B.chr[a])}:${mb(B.start[a])} – ${chrLabel(B.chr[b])}:${mb(B.end[b])}`;
}
