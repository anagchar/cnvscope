import { forwardRef, memo, useEffect, useImperativeHandle, useMemo, useRef, useState } from "react";
import { fetchCell, type CellDetail, type Meta } from "../api/client";
import { FONT, axes, chrLabel, css, esc, fit, fmt, fmtBp, fmtCount, fmtInt, niceTicks, tip } from "../lib/canvas";
import { binAt, binRange, chromLabels } from "../lib/genome";
import { fails } from "../lib/filters";
import { fmtMetric, isOutlier } from "../lib/metrics";
import { PALETTES } from "../lib/palettes";
import { discreteColor, isCall, scaleFor } from "../lib/scales";
import { useDerived } from "../state/derived";
import { ALLELE, useStore } from "../state/store";
import { FigureDialog } from "./FigureDialog";
import { DlIcon, useSize, useWidth } from "./widgets";

const PR = { l: 50, t: 12, r: 10, b: 26 };
type Vals = (number | null)[];

// recent cells, shared by every card that shows the selected cell
const cellCache = new WeakMap<Meta, Map<number, Promise<CellDetail>>>();
function cachedCell(meta: Meta, index: number): Promise<CellDetail> {
  let m = cellCache.get(meta);
  if (!m) cellCache.set(meta, (m = new Map()));
  let p = m.get(index);
  if (!p) {
    p = fetchCell(index);
    p.catch(() => m.delete(index));
    m.set(index, p);
    if (m.size > 32) m.delete(m.keys().next().value!);
  }
  return p;
}

/** The selected cell's layers, fetched when the selection changes. */
export function useCellDetail(index: number | null): CellDetail | null {
  const { meta } = useDerived();
  const [cell, setCell] = useState<{ i: number; c: CellDetail } | null>(null);
  useEffect(() => {
    if (index == null) return;
    let live = true;
    cachedCell(meta, index).then((c) => live && setCell({ i: index, c })).catch(() => live && setCell(null));
    return () => { live = false; };
  }, [index, meta]);
  return index != null && cell?.i === index ? cell.c : null;
}

function drawChromosomes(g: CanvasRenderingContext2D, meta: Meta, xS: (v: number) => number, gx0: number, gx1: number, top: number, h: number, labelY: number | null) {
  g.textAlign = "center"; g.textBaseline = "alphabetic"; g.font = `11px ${FONT}`;
  const spans: [string, number, number][] = [];
  for (const { chr, start_cum, end_cum } of meta.chr_info) {
    if (end_cum < gx0 || start_cum > gx1) continue;
    const a = xS(Math.max(start_cum, gx0)), b = xS(Math.min(end_cum, gx1));
    if (start_cum > gx0) { g.strokeStyle = css("--line2"); g.lineWidth = 1; g.beginPath(); g.moveTo(Math.round(a) + 0.5, top); g.lineTo(Math.round(a) + 0.5, top + h); g.stroke(); }
    spans.push([chr, a, b]);
  }
  if (labelY == null) return;
  g.fillStyle = css("--ink2");
  for (const [t, x] of chromLabels(g, spans)) g.fillText(t, x, labelY);
}

/** Runs of one state within a chromosome: [first bin, last bin + 1, state]. */
function runsOf(chr: string[], vals: Vals): [number, number, number][] {
  const out: [number, number, number][] = [], nb = vals.length;
  for (let s = 0; s < nb; ) {
    let e = s + 1;
    while (e < nb && vals[e] === vals[s] && chr[e] === chr[s]) e++;
    if (vals[s] != null) out.push([s, e, vals[s] as number]);
    s = e;
  }
  return out;
}

/* ── profile: copy-number calls over the bins, following the heatmap's genomic zoom ── */
const ProfilePlot = memo(forwardRef<HTMLCanvasElement, { cell: CellDetail; height?: number }>(function ProfilePlot({ cell, height = 300 }, ref) {
  const d = useDerived();
  const { meta, totalLayer, binsLayer, sem } = d;
  const view = useStore((s) => s.view);
  const zoom = useStore((s) => s.zoom);
  const showCi = useStore((s) => s.showCi);
  const theme = useStore((s) => s.theme);
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  useImperativeHandle(ref, () => cv.current!, []);
  const width = useWidth(box);
  const B = meta.bins, nb = meta.n_bins, genomeMax = B.end_cum[nb - 1];
  const [gx0, gx1] = zoom?.x ?? [0, genomeMax];

  const L = cell.layers;
  const isAS = view === ALLELE && !!L.major_cn && !!L.minor_cn;
  // calls: the integer total, or the view's own call layer when there is none
  const callLayer = totalLayer ?? (isCall(sem[view] ?? "count") ? view : null);
  const prof = callLayer ? L[callLayer] : null;
  const raw = useMemo(() => {
    if (!binsLayer || !L[binsLayer[0]]) return !callLayer && L[view] ? L[view] : null;
    const v = L[binsLayer[0]];
    if (binsLayer[1] === "cn" || !prof) return v;
    let s = 0, n = 0;
    for (const x of prof) if (x != null) { s += x; n++; }
    const psi = n ? s / n : 2;
    return v.map((x) => (x == null ? null : psi * Math.pow(2, x)));
  }, [L, binsLayer, prof, callLayer, view]);
  const ivMaj = meta.intervals.major_cn, ivMin = meta.intervals.minor_cn;
  const hasCi = isAS && showCi && !!ivMaj && !!ivMin && L[ivMaj[0]]?.some((v) => v != null);
  const cnScale = !!callLayer && (sem[callLayer] === "integer_cn" || sem[callLayer] === "allele_cn");
  const ternary = !!callLayer && sem[callLayer] === "ternary_state";
  const palette = useStore((s) => s.palette);
  const stateScale = useMemo(() => scaleFor("integer_cn", { palette: PALETTES[palette], dark: theme === "dark" }), [palette, theme]);
  const ternaryScale = useMemo(() => scaleFor("ternary_state", { palette: PALETTES[palette], dark: theme === "dark" }), [palette, theme]);

  useEffect(() => {
    const c = cv.current;
    if (!c) return;
    const g = fit(c, width, height);
    g.clearRect(0, 0, width, height);
    const plotW = width - PR.l - PR.r, plotH = height - PR.t - PR.b;
    const xS = (v: number) => PR.l + ((v - gx0) / (gx1 - gx0)) * plotW;
    const [b0, b1] = binRange(B.start_cum, B.end_cum, gx0, gx1);
    const nz = b1 - b0;

    let Y: (v: number) => number, ymax = 0;
    if (cnScale) {
      const top = (arr: Vals) => { let m = 0; for (const v of arr) if (v != null && v > m) m = v; return m; };
      ymax = Math.max(5, Math.min(12, (isAS ? top(L.major_cn) : top(prof!)) + 1));
      Y = (v) => PR.t + plotH - (Math.max(-0.3, Math.min(ymax + 0.4, v)) / (ymax + 0.4)) * plotH;
      g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
      for (let y = 0; y <= ymax; y++) {
        g.setLineDash(y === 2 ? [] : [3, 4]); g.beginPath(); g.moveTo(PR.l, Y(y)); g.lineTo(PR.l + plotW, Y(y)); g.stroke();
        g.fillText(String(y), PR.l - 8, Y(y));
      }
      g.setLineDash([]);
    } else if (ternary) {
      // gain / normal / loss calls (epiAneufinder)
      Y = (v) => PR.t + plotH - ((Math.max(-1.5, Math.min(1.5, v)) + 1.5) / 3) * plotH;
      g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
      ["Loss", "Normal", "Gain"].forEach((t, k) => {
        g.setLineDash(k === 1 ? [] : [3, 4]); g.beginPath(); g.moveTo(PR.l, Y(k - 1)); g.lineTo(PR.l + plotW, Y(k - 1)); g.stroke();
        g.fillText(t, PR.l - 8, Y(k - 1));
      });
      g.setLineDash([]);
    } else {
      // a generic layer: its own range
      const vals = (raw ?? prof ?? []).slice(b0, b1).filter((v): v is number => v != null).sort((a, b) => a - b);
      let lo = vals.length ? vals[Math.floor(0.005 * (vals.length - 1))] : 0, hi = vals.length ? vals[Math.floor(0.995 * (vals.length - 1))] : 1;
      if (prof) for (let i = b0; i < b1; i++) { const v = prof[i]; if (v != null) { lo = Math.min(lo, v); hi = Math.max(hi, v); } }
      if (!(hi > lo)) { lo -= 1; hi += 1; }
      const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;
      Y = (v) => PR.t + plotH - ((Math.max(lo, Math.min(hi, v)) - lo) / (hi - lo)) * plotH;
      g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
      for (const t of niceTicks(lo, hi, 5)) { g.beginPath(); g.moveTo(PR.l, Y(t)); g.lineTo(PR.l + plotW, Y(t)); g.stroke(); g.fillText(String(t), PR.l - 8, Y(t)); }
    }
    drawChromosomes(g, meta, xS, gx0, gx1, PR.t, plotH, PR.t + plotH + 18);
    g.save(); g.beginPath(); g.rect(PR.l, 0, plotW, height); g.clip();
    if (raw) {
      // bins coloured by their copy-number state, as in scgenome / SIGNALS profiles
      const pointCol = css("--points");
      g.globalAlpha = 0.8;
      const r = nz > 1500 ? 1.4 : nz > 400 ? 1.9 : 2.7;
      for (let i = b0; i < b1; i++) {
        const v = raw[i];
        if (v == null || !Number.isFinite(v)) continue;
        g.fillStyle = (cnScale && prof && discreteColor(stateScale, prof[i] ?? NaN)) || pointCol;
        g.beginPath(); g.arc(xS((B.start_cum[i] + B.end_cum[i]) / 2), Y(v), r, 0, 7); g.fill();
      }
      g.globalAlpha = 1;
    }
    const cap = (v: number) => (cnScale ? Math.min(v, ymax) : v);
    const runs = prof ? runsOf(B.chr, prof) : [];
    if (prof && ternary) {
      g.lineWidth = 4; g.lineCap = "round";
      for (const [s, e, v] of runs) {
        g.strokeStyle = discreteColor(ternaryScale, v) || css("--call");
        g.beginPath(); g.moveTo(xS(B.start_cum[s]) + 1, Y(v)); g.lineTo(xS(B.end_cum[e - 1]) - 1, Y(v)); g.stroke();
      }
    } else if (prof && !isAS) {
      g.strokeStyle = css("--call"); g.lineWidth = 3; g.lineCap = "round"; g.beginPath();
      for (const [s, e, v] of runs) { g.moveTo(xS(B.start_cum[s]) + 1, Y(cap(v))); g.lineTo(xS(B.end_cum[e - 1]) - 1, Y(cap(v))); }
      g.stroke();
    } else if (isAS) {
      const sets: [Vals, string, [string, string] | undefined, number][] = [[L.major_cn, "--major", ivMaj, 0.07], [L.minor_cn, "--minor", ivMin, -0.07]];
      for (const [arr, col, iv, off] of sets) {
        if (hasCi && iv) {
          const lo = L[iv[0]], hi = L[iv[1]];
          g.fillStyle = css(col); g.globalAlpha = 0.2;
          for (let i = b0; i < b1; i++) {
            const a = arr[i], l = lo[i], h = hi[i];
            if (a == null || l == null || h == null) continue;
            const yt = Y(Math.max(a, h)), yb = Y(Math.min(a, l));
            g.fillRect(xS(B.start_cum[i]), yt, xS(B.end_cum[i]) - xS(B.start_cum[i]), yb - yt);
          }
          g.globalAlpha = 1;
        }
        g.strokeStyle = css(col); g.lineWidth = 3; g.lineCap = "round"; g.beginPath();
        for (const [s, e, v] of runsOf(B.chr, arr)) { g.moveTo(xS(B.start_cum[s]) + 1, Y(Math.min(v, ymax) + off)); g.lineTo(xS(B.end_cum[e - 1]) - 1, Y(Math.min(v, ymax) + off)); }
        g.stroke();
      }
    }
    if (raw && prof) {
      g.strokeStyle = css("--mean"); g.lineWidth = 1.2; g.lineCap = "butt"; g.beginPath();
      for (const [s, e] of runs) {
        let m = 0, k = 0;
        for (let i = s; i < e; i++) { const v = raw[i]; if (v != null && Number.isFinite(v)) { m += v; k++; } }
        if (!k) continue;
        m /= k;
        g.moveTo(xS(B.start_cum[s]), Y(m)); g.lineTo(xS(B.end_cum[e - 1]), Y(m));
      }
      g.stroke();
    }
    g.restore();
    if (!raw && prof && !ternary) { g.fillStyle = css("--ink3"); g.textAlign = "left"; g.textBaseline = "alphabetic"; g.font = `11px ${FONT}`; g.fillText("no per-bin values in this dataset; showing the calls", PR.l + 8, PR.t + 14); }
  }, [cell, width, height, gx0, gx1, isAS, raw, prof, hasCi, cnScale, L, meta, B, nb, ivMaj, ivMin, theme, stateScale, ternary, ternaryScale]);

  const move = (e: React.MouseEvent) => {
    const r = cv.current!.getBoundingClientRect(), x = e.clientX - r.left, plotW = width - PR.l - PR.r;
    if (x < PR.l || x > PR.l + plotW) return tip(null);
    const b = binAt(B.end_cum, gx0 + ((x - PR.l) / plotW) * (gx1 - gx0));
    const line = isAS
      ? `major|minor <b>${L.major_cn[b] ?? "–"}|${L.minor_cn[b] ?? "–"}</b>`
      : prof ? `${esc(callLayer === totalLayer ? "copy number" : callLayer!)} <b>${prof[b] ?? "–"}</b>` : "";
    const rd = d.countLayer ? L[d.countLayer]?.[b] : null;
    tip(e, `<span class="k">${chrLabel(B.chr[b])}:${fmtBp(B.start[b])}–${fmtBp(B.end[b])}</span>${line ? "<br>" + line : ""}` +
      (raw && raw[b] != null ? `<br><span class="k">${cnScale ? "bin on the CN scale" : esc(view)}</span> ${fmt(raw[b])}` : "") +
      (rd != null ? `<br><span class="k">reads</span> ${rd}` : ""));
  };

  const items: [string, string, string][] = [];
  if (raw) items.push(["dot", cnScale ? "--points-key" : "--points", cnScale ? "bin, on the copy-number scale, coloured by its state" : `bin (${view})`]);
  if (raw && prof) items.push(["line", "--mean", "segment mean"]);
  if (isAS) items.push(["thick", "--major", "major allele"], ["thick", "--minor", "minor allele"]);
  else if (prof) items.push(["thick", "--call", callLayer === totalLayer ? "integer copy number" : callLayer!]);
  return (
    <div ref={box} style={{ width: "100%" }}>
      <canvas ref={cv} onMouseMove={move} onMouseLeave={() => tip(null)} />
      <div className="legend lines">
        {items.map(([k, c, t]) => (
          <span className="sw" key={t}><i style={{ width: k === "dot" ? 8 : 20, height: k === "thick" ? 4 : k === "dot" ? 8 : 2, borderRadius: k === "dot" ? 4 : 2, background: `var(${c})` }} />{t}</span>
        ))}
        {hasCi && <span className="sw"><i style={{ width: 20, height: 8, borderRadius: 2, background: "var(--major)", opacity: 0.25 }} />BAF 5–95% interval</span>}
      </div>
    </div>
  );
}));

/* ── read depth: per-bin read counts on the profile's x scale ── */
function ReadDepthPlot({ reads, height = 120 }: { reads: Vals; height?: number }) {
  const { meta } = useDerived();
  const zoom = useStore((s) => s.zoom);
  const theme = useStore((s) => s.theme);
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const width = useWidth(box);
  const B = meta.bins, genomeMax = B.end_cum[meta.n_bins - 1];
  const [gx0, gx1] = zoom?.x ?? [0, genomeMax];
  const T = 18, BOT = 6;
  useEffect(() => {
    const c = cv.current;
    if (!c) return;
    const g = fit(c, width, height);
    g.clearRect(0, 0, width, height);
    const plotW = width - PR.l - PR.r, plotH = height - T - BOT;
    const xS = (v: number) => PR.l + ((v - gx0) / (gx1 - gx0)) * plotW;
    const [b0, b1] = binRange(B.start_cum, B.end_cum, gx0, gx1);
    // the axis stops at the 99.5th percentile so a few extreme bins don't flatten the rest;
    // clipped bars get a mark on top
    let sum = 0, n = 0;
    const seen: number[] = [];
    for (let i = b0; i < b1; i++) { const v = reads[i]; if (v == null) continue; seen.push(v); sum += v; n++; }
    seen.sort((a, b) => a - b);
    const max = Math.max(1, seen.length ? seen[Math.floor(0.995 * (seen.length - 1))] * 1.1 : 1);
    const Y = (v: number) => T + plotH - (Math.min(v, max) / max) * plotH;
    g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
    for (const t of niceTicks(0, max, 3)) { g.beginPath(); g.moveTo(PR.l, Y(t)); g.lineTo(PR.l + plotW, Y(t)); g.stroke(); g.fillText(String(t), PR.l - 8, Y(t)); }
    drawChromosomes(g, meta, xS, gx0, gx1, T, plotH, null);
    g.save(); g.beginPath(); g.rect(PR.l, 0, plotW, height); g.clip();
    g.fillStyle = css("--series1"); g.globalAlpha = 0.75;
    for (let i = b0; i < b1; i++) {
      const v = reads[i];
      if (v == null) continue;
      const x0 = xS(B.start_cum[i]), x1 = xS(B.end_cum[i]), y = Y(v);
      g.fillRect(x0, y, Math.max(0.6, x1 - x0), T + plotH - y);
    }
    g.globalAlpha = 1;
    g.fillStyle = css("--serious");
    for (let i = b0; i < b1; i++) {
      const v = reads[i];
      if (v != null && v > max) g.fillRect(xS(B.start_cum[i]), T - 2, Math.max(1.5, xS(B.end_cum[i]) - xS(B.start_cum[i])), 2);
    }
    if (n) { g.strokeStyle = css("--ink3"); g.setLineDash([4, 4]); g.beginPath(); g.moveTo(PR.l, Y(sum / n)); g.lineTo(PR.l + plotW, Y(sum / n)); g.stroke(); g.setLineDash([]); }
    g.restore();
    g.fillStyle = css("--ink2"); g.textAlign = "left"; g.textBaseline = "alphabetic";
    g.fillText("Reads per bin", PR.l, 12);
    if (n) { g.fillStyle = css("--ink3"); g.textAlign = "right"; g.fillText(`dashed: mean ${fmt(sum / n, 1)}`, PR.l + plotW, 12); }
  }, [reads, width, height, gx0, gx1, meta, B, theme]);
  return <div ref={box} style={{ width: "100%", marginTop: 10 }}><canvas ref={cv} /></div>;
}

/* ── the card: status, metric pills, profile and read depth ── */
export function ProfileCard() {
  const d = useDerived();
  const { meta, vals, metAll, passing, reasons, countLayer } = d;
  const selected = useStore((s) => s.selected);
  const view = useStore((s) => s.view);
  const cell = useCellDetail(selected);
  const plot = useRef<HTMLCanvasElement>(null);
  const [fig, setFig] = useState(false);
  const name = selected != null ? meta.cells[selected] : null;
  const valOf = (k: string) => (selected != null && vals[k] ? vals[k][selected] : NaN);
  const failsFilter = (metric: string) => fails(d.filters[metric], valOf(metric));
  const metricOf = (k: string) => metAll.find((m) => m.key === k);

  const pills: [string, string, string, string?][] = [];
  if (selected != null) {
    const add = (k: string, label: string, v: string, title?: string) => { if (!Number.isNaN(valOf(k))) pills.push([k, label, v, title]); };
    add("ploidy", "ploidy", fmt(valOf("ploidy")));
    add("n_reads_total", "reads", fmtCount(valOf("n_reads_total")), `${fmtInt(valOf("n_reads_total"))} reads in this cell's bins`);
    add("reads_per_bin", "reads / bin", fmt(valOf("reads_per_bin"), 1));
    add("coverage", "coverage", fmtMetric("coverage", valOf("coverage")));
    add("mapd", "MAPD", fmt(valOf("mapd"), 3));
    add("median_residual", "residual", fmt(valOf("median_residual"), 3));
    add("gini", "Gini", fmt(valOf("gini"), 3));
    add("bins_with_cna", "CNA bins", fmtInt(valOf("bins_with_cna")), `bins whose copy number is not 2, of ${meta.n_bins}`);
    add("breakpoints", "breakpoints", fmtInt(valOf("breakpoints")));
  }
  const ct = selected != null ? meta.obs.cell_type?.[selected] : null;
  const status = selected == null ? null : passing[selected]
    ? <span className="status" style={{ color: "var(--good)" }}><svg width="12" height="12" viewBox="0 0 12 12"><circle cx="6" cy="6" r="6" fill="currentColor" /><path d="M3.4 6.2 5.2 8l3.4-3.6" stroke="#fff" strokeWidth="1.5" fill="none" /></svg>passes QC</span>
    : <span className="status" style={{ color: "var(--serious)" }}><svg width="12" height="12" viewBox="0 0 12 12"><path d="M6 .8 11.4 10.6H.6z" fill="currentColor" /><path d="M6 4.3v3M6 8.8v.2" stroke="#fff" strokeWidth="1.4" /></svg>flagged{reasons[selected] ? ": " + reasons[selected] : ""}</span>;
  const reads = cell && countLayer ? cell.layers[countLayer] : null;

  return (
    <div className="card">
      {fig && <FigureDialog kind="profile" onClose={() => setFig(false)} />}
      <div className="hd">
        <h2>{name ? <span className="mono">{name}</span> : "Cell profile"}</h2>
        <span className="sub">{status}{ct && <span className="ok"> · {String(ct)}</span>}{view === ALLELE && d.hasAS && <span className="ok"> · allele-specific</span>}</span>
        <button className="iconbtn" onClick={() => setFig(true)} disabled={selected == null} title="Download a publication figure"><DlIcon />Figure</button>
      </div>
      {pills.length > 0 && (
        <div className="pills">
          {pills.map(([k, label, v, title]) => (
            <div className={"pill" + (failsFilter(k) ? " bad" : "")} key={k} title={title || metricOf(k)?.help}><span className="k">{label}</span><span className="v">{v}</span></div>
          ))}
        </div>
      )}
      <div className="bd">
        {selected == null ? <div className="empty">Click a row of the heatmap, or pick a cell on the left.</div>
          : !cell ? <div className="empty">Loading…</div>
          : <>
            <ProfilePlot ref={plot} cell={cell} />
            {reads && <ReadDepthPlot reads={reads} />}
          </>}
      </div>
    </div>
  );
}

/* ── Lorenz curve of the reads over bins ── */
export function LorenzCard() {
  const { meta, vals, countLayer } = useDerived();
  const selected = useStore((s) => s.selected);
  const theme = useStore((s) => s.theme);
  const cell = useCellDetail(selected);
  const reads = cell && countLayer ? cell.layers[countLayer] : null;
  const gini = selected != null ? vals.gini?.[selected] : undefined;
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  const { w: width, h: height } = useSize(box, 260);
  const points = useMemo(() => {
    if (!reads) return null;
    const v = reads.filter((x): x is number => x != null && x >= 0).sort((a, b) => a - b);
    const total = v.reduce((a, b) => a + b, 0);
    if (!v.length || !total) return null;
    const pts: [number, number][] = [[0, 0]];
    let cum = 0;
    v.forEach((x, i) => { cum += x; pts.push([(i + 1) / v.length, cum / total]); });
    return pts;
  }, [reads]);
  useEffect(() => {
    const c = cv.current;
    if (!c || !points) return;
    const g = fit(c, width, height);
    g.clearRect(0, 0, width, height);
    const b = { l: 52, t: 14, w: width - 66, h: height - 58 };
    const ticks = [0, 0.25, 0.5, 0.75, 1], pct = (v: number) => Math.round(v * 100) + "%";
    const { tx, ty } = axes(g, b, ticks, ticks, { x0: 0, x1: 1, y0: 0, y1: 1, xlab: "share of bins, fewest reads first", ylab: "share of reads", xfmt: pct, yfmt: pct });
    g.strokeStyle = css("--ink3"); g.setLineDash([4, 4]); g.beginPath(); g.moveTo(tx(0), ty(0)); g.lineTo(tx(1), ty(1)); g.stroke(); g.setLineDash([]);
    g.fillStyle = css("--accent-soft"); g.beginPath(); g.moveTo(tx(0), ty(0));
    for (const [x, y] of points) g.lineTo(tx(x), ty(y));
    g.lineTo(tx(1), ty(1)); g.closePath(); g.fill();
    g.strokeStyle = css("--accent"); g.lineWidth = 2; g.lineJoin = "round"; g.beginPath();
    points.forEach(([x, y], i) => (i ? g.lineTo(tx(x), ty(y)) : g.moveTo(tx(x), ty(y))));
    g.stroke();
    if (gini != null && !Number.isNaN(gini)) { g.fillStyle = css("--ink"); g.font = `600 12px ${FONT}`; g.textAlign = "left"; g.textBaseline = "top"; g.fillText(`Gini ${gini.toFixed(3)}`, b.l + 10, b.t + 6); }
  }, [points, width, height, gini, theme]);
  void meta;
  return (
    <div className="card">
      <div className="hd"><h2>Coverage evenness</h2><span className="sub">Lorenz curve of the reads over bins</span></div>
      <div className="bd">
        <div ref={box} className="fill" style={{ minHeight: points ? 260 : 0 }}>
          {points ? <canvas ref={cv} /> : <div className="empty">No per-bin read counts in this dataset.</div>}
        </div>
      </div>
    </div>
  );
}

/* ── this cell's QC against all cells ── */
export function CellQCCard() {
  const { metAll, vals, filters } = useDerived();
  const selected = useStore((s) => s.selected);
  const valOf = (k: string) => (selected != null && vals[k] ? vals[k][selected] : NaN);
  return (
    <div className="card">
      <div className="hd"><h2>This cell's QC</h2><span className="sub">against all cells · ▲ beyond median ± 3 MAD on the bad side</span></div>
      <div className="bd tablewrap">
        {metAll.length ? (
          <table>
            <thead><tr><th>metric</th><th>cell</th><th>median</th><th /></tr></thead>
            <tbody>
              {metAll.map((m) => {
                const v = valOf(m.key), out = isOutlier(m, v), fail = fails(filters[m.key], v);
                return (
                  <tr key={m.key}>
                    <td title={m.help}>{m.label}</td>
                    <td style={{ fontWeight: 600, color: fail ? "var(--serious)" : undefined }}>{fmtMetric(m.key, v)}</td>
                    <td className="ok">{fmtMetric(m.key, m.median)}</td>
                    <td style={{ color: "var(--serious)", width: 1, whiteSpace: "nowrap" }}>{fail ? "filter" : out ? "▲" : ""}</td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        ) : <div className="empty">This dataset has no per-cell QC metrics.</div>}
      </div>
    </div>
  );
}
