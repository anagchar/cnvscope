import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { HeatmapGL } from "../gl/heatmapGL";
import { CATEGORICAL, SEQUENTIAL } from "../lib/palettes";
import { FONT, chrLabel, css, esc, fit, fmt, fmtBp, rgbOf, swatch, tip } from "../lib/canvas";
import { layoutVisible } from "../lib/dendrogram";
import { STAIN_COLORS, binAt, binColumns, chromLabels, useCytobands } from "../lib/genome";
import { buildLut, type Scale } from "../lib/scales";
import { useDerived } from "../state/derived";
import { ALLELE, useStore, type Zoom } from "../state/store";
import { FigureDialog } from "./FigureDialog";
import { DlIcon, useWidth } from "./widgets";

// Rows of the heatmap start below the chromosome labels and the ideogram
const HM_TOP = 40, HM_BOT = 4, LABEL_Y = 13, IDEO_Y = 22, IDEO_H = 12;
const DENDRO_W = 92, ANNOT_W = 12, QC_W = 6;

export function useCategoryColors(values: (string | number | boolean | null)[] | undefined) {
  return useMemo(() => {
    const map = new Map<string, string>();
    if (!values) return map;
    const keys = [...new Set(values.filter((v) => v != null).map(String))].sort();
    keys.forEach((k, i) => map.set(k, CATEGORICAL[i % CATEGORICAL.length]));
    return map;
  }, [values]);
}

export function HeatmapCard() {
  const d = useDerived();
  const { meta, shown, rowOf, heat, dendro } = d;
  const view = useStore((s) => s.view);
  const zoom = useStore((s) => s.zoom);
  const set = useStore((s) => s.set);
  const selected = useStore((s) => s.selected);
  const theme = useStore((s) => s.theme);
  const height = useStore((s) => s.heatmapH);
  const showDendro = useStore((s) => s.showDendro) && !!dendro;
  const annotation = useStore((s) => s.annotation);
  const filtersOpen = useStore((s) => s.filtersOpen);
  const nFlagged = meta.n_cells - d.nPass;
  const qcRef = useRef<HTMLCanvasElement>(null);
  const cytobands = useCytobands(meta.genome);

  const nShown = shown.length;
  const nb = meta.n_bins;
  const B = meta.bins;
  const genomeMax = nb ? B.end_cum[nb - 1] : 1;
  const plotH = height - HM_TOP - HM_BOT;
  const [gx0, gx1] = zoom?.x ?? [0, genomeMax];
  const [y0, y1] = zoom?.y ?? [0, nShown];
  const xSpan = gx1 - gx0 || 1, ySpan = y1 - y0 || 1;

  const boxRef = useRef<HTMLDivElement>(null);
  const cw = useWidth(boxRef);
  const glRef = useRef<HTMLCanvasElement>(null);
  const frameRef = useRef<HTMLCanvasElement>(null);
  const overRef = useRef<HTMLCanvasElement>(null);
  const dendroRef = useRef<HTMLCanvasElement>(null);
  const annotRef = useRef<HTMLCanvasElement>(null);
  const gl = useRef<HeatmapGL | null>(null);
  const [glError, setGlError] = useState<string | null>(null);
  const [fig, setFig] = useState(false);

  const toPxX = useCallback((gx: number) => ((gx - gx0) / xSpan) * cw, [gx0, xSpan, cw]);
  const toPxY = useCallback((ci: number) => HM_TOP + ((ci - y0) / ySpan) * plotH, [y0, ySpan, plotH]);
  const toGx = (px: number) => gx0 + (px / cw) * xSpan;
  const toCi = (py: number) => y0 + ((py - HM_TOP) / plotH) * ySpan;

  // ── GPU layer ──
  useEffect(() => {
    if (!glRef.current || gl.current) return;
    try {
      gl.current = new HeatmapGL(glRef.current);
    } catch (e) {
      setGlError(String(e instanceof Error ? e.message : e));
    }
  }, []);
  useEffect(() => {
    const r = gl.current;
    if (!r || !heat.matrix) return;
    try {
      r.setData(heat.matrix, meta.n_cells, nb);
      setGlError(null);
    } catch (e) {
      setGlError(String(e instanceof Error ? e.message : e));
    }
  }, [heat.matrix, meta, nb]);
  useEffect(() => { gl.current?.setRows(shown); }, [shown]);
  useEffect(() => {
    const r = gl.current, cv = glRef.current;
    if (!r || !cv || !heat.matrix || !nShown) return;
    const dpr = Math.max(1, devicePixelRatio || 1);
    const W = Math.max(1, Math.round(cw * dpr)), H = Math.max(1, Math.round(plotH * dpr));
    cv.width = W; cv.height = H;
    cv.style.width = cw + "px"; cv.style.height = plotH + "px";
    r.setCols(binColumns(B, gx0, gx0 + xSpan, W) as Float32Array);
    r.setScale(heat.scale, buildLut(heat.scale, rgbOf));
    r.draw(y0, y1, rgbOf("--cnna"), rgbOf("--panel"));
  }, [heat, cw, plotH, gx0, xSpan, y0, y1, nShown, B, theme, shown]);

  // ── frame: labels, ideogram, chromosome separators ──
  useEffect(() => {
    const cv = frameRef.current;
    if (!cv) return;
    const g = fit(cv, cw, height);
    g.clearRect(0, 0, cw, height);
    if (!nShown) {
      g.fillStyle = css("--ink3"); g.font = `13px ${FONT}`; g.fillText("No cell passes the filters.", 12, HM_TOP + 30);
      return;
    }
    const dark = document.documentElement.dataset.theme === "dark";
    if (cytobands?.length) {
      const off: Record<string, number> = {};
      meta.chr_info.forEach((c) => (off[c.chr] = c.start_cum));
      g.fillStyle = dark ? css("--field") : "#ffffff";
      g.fillRect(0, IDEO_Y, cw, IDEO_H);
      for (const band of cytobands) {
        const o = off[band.chr];
        if (o === undefined) continue;
        const s = o + band.start, e = o + band.end;
        if (e <= gx0 || s >= gx1) continue;
        const a = Math.max(0, Math.floor(toPxX(Math.max(s, gx0)))), b = Math.min(cw, Math.ceil(toPxX(Math.min(e, gx1))));
        if (b <= a) continue;
        g.fillStyle = band.stain === "gneg" && dark ? css("--field") : STAIN_COLORS[band.stain] || "#AAAAAA";
        g.fillRect(a, IDEO_Y, b - a, IDEO_H);
      }
      g.strokeStyle = css("--line2");
      g.lineWidth = 0.5;
      for (const { start_cum, end_cum } of meta.chr_info) {
        if (end_cum <= gx0 || start_cum >= gx1) continue;
        const a = Math.max(0, Math.floor(toPxX(Math.max(start_cum, gx0)))), b = Math.min(cw, Math.ceil(toPxX(Math.min(end_cum, gx1))));
        if (b > a + 1) g.strokeRect(a + 0.5, IDEO_Y + 0.5, b - a - 1, IDEO_H - 1);
      }
    }
    g.fillStyle = css("--hm-gap");
    for (const { start_cum } of meta.chr_info) {
      if (start_cum <= gx0 || start_cum > gx1) continue;
      g.fillRect(toPxX(start_cum) - 0.75, HM_TOP, 1.5, plotH);
    }
    g.fillStyle = css("--ink2"); g.font = `11px ${FONT}`; g.textAlign = "center"; g.textBaseline = "alphabetic";
    const spans = meta.chr_info
      .filter(({ start_cum, end_cum }) => end_cum >= gx0 && start_cum <= gx1)
      .map(({ chr, start_cum, end_cum }) => [chr, toPxX(Math.max(start_cum, gx0)), toPxX(Math.min(end_cum, gx1))] as [string, number, number]);
    for (const [t, x] of chromLabels(g, spans)) g.fillText(t, x, LABEL_Y);
  }, [cw, height, plotH, meta, gx0, gx1, toPxX, cytobands, theme, nShown]);

  // ── overlay: selected cell, brush ──
  const brush = useRef<{ sx: number; sy: number; cx: number; cy: number } | null>(null);
  const dragging = useRef(false);
  const drawOverlay = useCallback(() => {
    const cv = overRef.current;
    if (!cv) return;
    const g = fit(cv, cw, height);
    g.clearRect(0, 0, cw, height);
    if (selected != null && rowOf[selected] >= 0) {
      const i = rowOf[selected], a = toPxY(i), b = toPxY(i + 1);
      if (b > HM_TOP && a < HM_TOP + plotH) {
        const hh = Math.max(2, b - a);
        g.strokeStyle = css("--accent"); g.lineWidth = 2;
        g.strokeRect(1, (a + b) / 2 - hh / 2 - 1, cw - 2, hh + 2);
      }
    }
    const br = brush.current;
    if (br && dragging.current) {
      const x = Math.min(br.sx, br.cx), y = Math.min(br.sy, br.cy), w = Math.abs(br.cx - br.sx), h = Math.abs(br.cy - br.sy);
      g.fillStyle = css("--accent-soft"); g.strokeStyle = css("--accent"); g.lineWidth = 1;
      g.fillRect(x, y, w, h); g.strokeRect(x, y, w, h);
    }
  }, [cw, height, selected, rowOf, toPxY, plotH]);
  useEffect(() => { drawOverlay(); }, [drawOverlay, theme]);

  // ── dendrogram ──
  const layout = useMemo(
    () => (dendro && showDendro && nShown > 1 ? layoutVisible(dendro.merge, dendro.height, rowOf, y0, y1, nShown) : null),
    [dendro, showDendro, rowOf, y0, y1, nShown],
  );
  useEffect(() => {
    const cv = dendroRef.current;
    if (!cv || !layout || !layout.maxH) return;
    const g = fit(cv, DENDRO_W, height), PAD = 6;
    g.clearRect(0, 0, DENDRO_W, height);
    const xs = (h: number) => DENDRO_W - PAD - Math.sqrt(h / layout.maxH) * (DENDRO_W - PAD * 2);
    g.strokeStyle = css("--ink3"); g.lineWidth = 0.8; g.beginPath();
    for (const s of layout.segments) {
      if (Math.max(s.y1, s.y2) < y0 || Math.min(s.y1, s.y2) > y1) continue;
      const a = s.y1 === s.y2 ? s.y1 : Math.max(y0, Math.min(y1, s.y1));
      const b = s.y1 === s.y2 ? s.y2 : Math.max(y0, Math.min(y1, s.y2));
      g.moveTo(xs(s.x1), toPxY(a)); g.lineTo(xs(s.x2), toPxY(b));
    }
    g.stroke();
  }, [layout, height, y0, y1, toPxY, theme]);

  // ── annotation track ──
  const annotValues = annotation ? meta.obs[annotation] : undefined;
  const annotKind = annotation ? meta.obs_kinds[annotation] : undefined;
  const catColors = useCategoryColors(annotKind === "categorical" ? annotValues : undefined);
  useEffect(() => {
    const cv = annotRef.current;
    if (!cv || !annotValues) return;
    const g = fit(cv, ANNOT_W, height);
    g.clearRect(0, 0, ANNOT_W, height);
    let lo = Infinity, hi = -Infinity;
    if (annotKind === "numeric") for (const v of annotValues) if (typeof v === "number") { lo = Math.min(lo, v); hi = Math.max(hi, v); }
    const accent = rgbOf(SEQUENTIAL[4]), base = rgbOf(SEQUENTIAL[0]);
    for (let r = Math.max(0, Math.floor(y0)); r < Math.min(nShown, Math.ceil(y1)); r++) {
      const v = annotValues[shown[r]];
      if (v == null) continue;
      if (annotKind === "numeric") {
        const t = hi > lo ? ((v as number) - lo) / (hi - lo) : 0.5;
        g.fillStyle = `rgb(${base.map((b, k) => Math.round(b + (accent[k] - b) * t)).join(",")})`;
      } else g.fillStyle = catColors.get(String(v)) ?? css("--cnna");
      const a = toPxY(r), b = toPxY(r + 1);
      g.fillRect(0, Math.max(HM_TOP, a), ANNOT_W, Math.max(1, b - a + 0.5));
    }
  }, [annotValues, annotKind, catColors, height, y0, y1, nShown, shown, toPxY, theme]);

  // ── QC strip: amber for each flagged cell shown ──
  useEffect(() => {
    const cv = qcRef.current;
    if (!cv) return;
    const g = fit(cv, QC_W, height);
    g.clearRect(0, 0, QC_W, height);
    g.fillStyle = css("--qc-other");
    g.fillRect(0, HM_TOP, QC_W, plotH);
    g.fillStyle = css("--qc-fail");
    for (let r = Math.max(0, Math.floor(y0)); r < Math.min(nShown, Math.ceil(y1)); r++) {
      if (d.passing[shown[r]]) continue;
      const a = toPxY(r), b = toPxY(r + 1);
      g.fillRect(0, Math.max(HM_TOP, a), QC_W, Math.max(1, b - a + 0.5));
    }
  }, [d.passing, shown, height, plotH, y0, y1, nShown, toPxY, theme]);

  // ── interactions ──
  const xy = (e: { clientX: number; clientY: number }) => {
    const r = overRef.current!.getBoundingClientRect();
    return { x: e.clientX - r.left, y: e.clientY - r.top };
  };
  const hover = (e: React.MouseEvent, p: { x: number; y: number }) => {
    const ci = Math.floor(toCi(p.y));
    if (p.y < HM_TOP || ci < 0 || ci >= nShown || !heat.matrix) return tip(null);
    const row = shown[ci], bi = binAt(B.end_cum, toGx(p.x));
    const v = heat.matrix[row * nb + bi];
    let val: string;
    if (view === ALLELE) {
      const k = Number.isNaN(v) ? null : heat.scale.labels[v];
      val = k ? `${swatch(heat.scale.stops[v])} ${k}` : "no value";
    } else if (Number.isNaN(v)) val = "no value";
    else if (heat.scale.kind === "discrete") {
      const k = Math.max(heat.scale.lo, Math.min(heat.scale.hi, Math.round(v))) - heat.scale.lo;
      val = `${swatch(css(heat.scale.stops[k]))} ${Number.isInteger(v) ? v : fmt(v, 2)}`;
    } else val = fmt(v, 3);
    const ct = meta.obs.cell_type?.[row];
    tip(e, `<b class="mono">${esc(meta.cells[row])}</b>${ct ? ` <span class="k">· ${esc(String(ct))}</span>` : ""}<br><span class="k">${chrLabel(B.chr[bi])}:${fmtBp(B.start[bi])}–${fmtBp(B.end[bi])}</span><br>${view === ALLELE ? "major|minor" : esc(heat.label)} <b>${val}</b>`);
  };
  const setZoom = (z: Zoom | null) => set({ zoom: z });
  const onDown = (e: React.MouseEvent) => {
    if (e.button !== 0) return;
    const p = xy(e);
    brush.current = { sx: p.x, sy: p.y, cx: p.x, cy: p.y };
    dragging.current = true;
    tip(null);
  };
  const onMove = (e: React.MouseEvent) => {
    const p = xy(e);
    if (!dragging.current) return hover(e, p);
    brush.current!.cx = p.x; brush.current!.cy = p.y;
    drawOverlay();
  };
  const onUp = (e: React.MouseEvent) => {
    if (!dragging.current) return;
    dragging.current = false;
    const p = xy(e), b = brush.current;
    brush.current = null;
    if (b) {
      const dx = Math.abs(b.sx - p.x), dy = Math.abs(b.sy - p.y);
      if (dx < 5 && dy < 5) {
        const ci = Math.floor(toCi(p.y));
        if (ci >= 0 && ci < nShown) set({ selected: shown[ci] });
      } else if (dx > 8 || dy > 8) {
        const a = toGx(Math.min(b.sx, p.x)), c = toGx(Math.max(b.sx, p.x));
        const r0 = toCi(Math.max(HM_TOP, Math.min(b.sy, p.y))), r1 = toCi(Math.max(b.sy, p.y));
        setZoom({ x: [Math.max(0, a), Math.min(genomeMax, c)], y: [Math.max(0, r0), Math.min(nShown, r1)] });
      }
    }
    drawOverlay();
  };
  // the wheel pans a zoomed heatmap (rows with the wheel, the genome sideways)
  useEffect(() => {
    const el = overRef.current;
    if (!el) return;
    const h = (e: WheelEvent) => {
      const yZ = y0 !== 0 || y1 !== nShown, xZ = gx0 !== 0 || gx1 !== genomeMax;
      if (!yZ && !xZ) return;
      e.preventDefault();
      let nx: [number, number] = [gx0, gx1], ny: [number, number] = [y0, y1];
      if (e.deltaY && yZ) {
        const span = y1 - y0, a = Math.max(0, Math.min(nShown - span, y0 + span * 0.1 * Math.sign(e.deltaY)));
        ny = [a, a + span];
      }
      if (e.deltaX && xZ) {
        const span = gx1 - gx0, a = Math.max(0, Math.min(genomeMax - span, gx0 + span * 0.1 * Math.sign(e.deltaX)));
        nx = [a, a + span];
      }
      set({ zoom: { x: nx, y: ny } });
    };
    el.addEventListener("wheel", h, { passive: false });
    return () => el.removeEventListener("wheel", h);
  }, [gx0, gx1, y0, y1, nShown, genomeMax, set]);

  const zoomLabel = () => {
    if (!zoom) return "";
    const a = binAt(B.end_cum, zoom.x[0]), b = binAt(B.end_cum, zoom.x[1]);
    return `${chrLabel(B.chr[a])}:${fmtBp(B.start[a])} – ${chrLabel(B.chr[b])}:${fmtBp(B.end[b])}`;
  };
  const msg = glError ?? (!heat.matrix ? `Loading ${heat.label}${heat.progress != null ? ` … ${Math.round(100 * heat.progress)}%` : "…"}` : null);

  return (
    <div className="card">
      {fig && <FigureDialog kind="heatmap" onClose={() => setFig(false)} />}
      <div className="hd">
        <h2>{meta.layers.find((l) => l.name === meta.default_layer)?.semantics === "ternary_state" ? "Genome-wide gains and losses" : "Genome-wide copy number"}</h2>
        <span className="sub">
          {nShown.toLocaleString()} of {meta.n_cells.toLocaleString()} cells · {nb.toLocaleString()} bins{zoom ? " (zoomed)" : ""}
          {view === ALLELE ? " · allele-specific" : view !== meta.default_layer ? ` · ${view}` : ""}
        </span>
        {zoom && <button className="zoomchip" onClick={() => setZoom(null)}>{zoomLabel()}  ✕</button>}
        <button className={"iconbtn" + (filtersOpen ? " active" : "")} onClick={() => set({ filtersOpen: !filtersOpen })} title="QC filters beside the heatmap" aria-pressed={filtersOpen}>
          <svg width="12" height="12" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.6" strokeLinecap="round"><path d="M2 4h12M4.5 8h7M7 12h2" /></svg>Filters{nFlagged ? ` · ${nFlagged.toLocaleString()} flagged` : ""}
        </button>
        <button className="iconbtn" onClick={() => setFig(true)} title="Download a publication figure"><DlIcon />Figure</button>
      </div>
      <div className="bd">
        <div className="hmwrap" style={{ display: "flex", width: "100%" }}>
          {showDendro && <canvas ref={dendroRef} style={{ width: DENDRO_W, height, flexShrink: 0 }} />}
          <canvas ref={qcRef} title="QC: amber = flagged by a filter" style={{ width: QC_W, height, flexShrink: 0, marginRight: 3, borderRadius: 3 }} />
          {annotValues && <canvas ref={annotRef} style={{ width: ANNOT_W, height, flexShrink: 0, marginRight: 3 }} />}
          <div ref={boxRef} style={{ flex: 1, minWidth: 0, position: "relative", height }}>
            <canvas ref={glRef} style={{ position: "absolute", left: 0, top: HM_TOP }} />
            <canvas ref={frameRef} style={{ position: "absolute", left: 0, top: 0 }} />
            <canvas
              ref={overRef}
              style={{ position: "absolute", left: 0, top: 0, cursor: "crosshair" }}
              onMouseDown={onDown}
              onMouseMove={onMove}
              onMouseUp={onUp}
              onMouseLeave={() => {
                tip(null);
                if (dragging.current) { dragging.current = false; brush.current = null; drawOverlay(); }
              }}
              onDoubleClick={() => setZoom(null)}
            />
            {msg && <div className="hm-msg" style={{ top: HM_TOP }}>{msg}</div>}
          </div>
        </div>
        <div className="keyrow">
          <span className="keytitle">{view === ALLELE ? "Major | minor copies" : heat.label === meta.default_layer ? (meta.layers.find((l) => l.name === meta.default_layer)?.semantics === "ternary_state" ? "Gain / loss" : "Total copy number") : heat.label}</span>
          <Legend scale={heat.scale} allele={view === ALLELE} />
          <span className="ok" style={{ marginLeft: "auto" }}>drag to zoom · double-click to reset</span>
        </div>
      </div>
    </div>
  );
}

const cssColor = (t: string) => (t.startsWith("--") ? `var(${t})` : t);

function Legend({ scale, allele }: { scale: Scale; allele: boolean }) {
  if (allele)
    return (
      <div className="legend as" aria-label="allele-specific state colours">
        {scale.labels.map((k, i) => <span className="sw" key={k}><i style={{ background: scale.stops[i] }} />{k}</span>)}
        <span className="sw"><i style={{ background: "var(--cnna)" }} />NA</span>
      </div>
    );
  if (scale.kind === "discrete")
    return (
      <div className="legend" aria-label="colours">
        {scale.stops.map((t, i) => <span className="sw" key={t + i}><i style={{ background: cssColor(t) }} />{scale.labels[i]}</span>)}
        <span className="sw"><i style={{ background: "var(--cnna)" }} />NA</span>
      </div>
    );
  return (
    <div className="legend ramp">
      <span>{scale.labels[0]}</span>
      <i style={{ background: `linear-gradient(90deg, ${scale.stops.map((t) => (t.startsWith("--") ? `var(${t})` : t)).join(",")})` }} />
      <span>{scale.labels[scale.labels.length - 1]}</span>
    </div>
  );
}
