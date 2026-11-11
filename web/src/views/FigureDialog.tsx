import { useEffect, useMemo, useRef, useState } from "react";
import { drawHeatmapFigure, drawProfileFigure, figTheme, regionLabel, type CommonOpts } from "../lib/figure";
import { useCytobands } from "../lib/genome";
import { CATEGORICAL, PALETTES } from "../lib/palettes";
import { alleleScale, scaleFor } from "../lib/scales";
import { useDerived } from "../state/derived";
import { ALLELE, useStore } from "../state/store";
import { useCellDetail } from "./ProfileCard";
import { Switch } from "./widgets";

export type FigureKind = "heatmap" | "profile";

const MM_PER_PT = 25.4 / 72;
const WIDTHS: [string, number][] = [["Single column · 89 mm", 89], ["1.5 column · 120 mm", 120], ["Double column · 183 mm", 183]];

function Row({ label, hint, children }: { label: string; hint?: string; children: React.ReactNode }) {
  return (
    <div className="frow">
      <div><label>{label}</label>{hint && <small>{hint}</small>}</div>
      <div className="fctl">{children}</div>
    </div>
  );
}

function Seg<T extends string | number>({ value, options, onChange }: { value: T; options: [T, string][]; onChange: (v: T) => void }) {
  return (
    <div className="seg">
      {options.map(([v, l]) => <button key={String(v)} className={v === value ? "on" : ""} onClick={() => onChange(v)}>{l}</button>)}
    </div>
  );
}

/** "Download figure": a live preview of a publication figure and its options. */
export function FigureDialog({ kind, onClose }: { kind: FigureKind; onClose: () => void }) {
  const d = useDerived();
  const { meta, heat, shown, rowOf, sem, totalLayer, binsLayer, countLayer, hasAS } = d;
  const st = useStore();
  const cell = useCellDetail(kind === "profile" ? st.selected : null);
  const cytobands = useCytobands(meta.genome);
  const genomeMax = meta.bins.end_cum[meta.n_bins - 1];
  const isAllele = st.view === ALLELE;
  const cellName = st.selected != null ? meta.cells[st.selected] : "";

  const [format, setFormat] = useState<"png" | "pdf">("png");
  const [widthMm, setWidthMm] = useState(183);
  const [heightMm, setHeightMm] = useState(kind === "heatmap" ? 120 : 75);
  const [dpi, setDpi] = useState(300);
  const [themeName, setThemeName] = useState<"print" | "screen">("print");
  const [fontPt, setFontPt] = useState(7);
  const [title, setTitle] = useState(kind === "heatmap" ? meta.name : cellName);
  const [showSubtitle, setShowSubtitle] = useState(true);
  const [legend, setLegend] = useState(true);
  const [region, setRegion] = useState<"genome" | "zoom">(st.zoom ? "zoom" : "genome");
  // heatmap
  const [dendro, setDendro] = useState(st.showDendro && !!d.dendro);
  const [ideogram, setIdeogram] = useState(true);
  const annotCols = Object.keys(meta.obs_kinds).filter((c) => meta.obs_kinds[c] === "categorical");
  const [annotation, setAnnotation] = useState<string>(st.annotation && annotCols.includes(st.annotation) ? st.annotation : annotCols.includes("cell_type") ? "cell_type" : "");
  // profile
  const [reads, setReads] = useState(!!countLayer);
  const [colorBins, setColorBins] = useState(true);
  const [ci, setCi] = useState(st.showCi);
  const [busy, setBusy] = useState(false);

  const zoomed = !!st.zoom;
  const gx: [number, number] = region === "zoom" && st.zoom ? st.zoom.x : [0, genomeMax];
  const rows: [number, number] = kind === "heatmap" && region === "zoom" && st.zoom ? st.zoom.y : [0, shown.length];
  const subtitle = !showSubtitle ? "" : kind === "heatmap"
    ? `${meta.tool} · ${(Math.round(rows[1] - rows[0])).toLocaleString()} cells · ${isAllele ? "allele-specific copy number" : "total copy number"}${region === "zoom" ? " · " + regionLabel(meta, gx) : ""}`
    : `${meta.tool}${st.selected != null && d.vals.ploidy ? ` · ploidy ${d.vals.ploidy[st.selected].toFixed(2)}` : ""}${isAllele && hasAS ? " · allele-specific" : ""}${region === "zoom" ? " · " + regionLabel(meta, gx) : ""}`;

  const theme = themeName === "print" ? figTheme(true) : figTheme(false);
  const dark = themeName === "screen" && st.theme === "dark";
  const palette = PALETTES[st.palette];
  // the scale again for the figure's theme (the neutral grey differs on dark backgrounds)
  const scale = useMemo(() => {
    if (isAllele) return alleleScale(palette, dark);
    const s = sem[st.view] ?? "integer_cn";
    return scaleFor(s, { palette, dark, countMax: s === "count" ? heat.scale.hi : undefined });
  }, [isAllele, palette, dark, sem, st.view, heat.scale.hi]);
  const stateScale = useMemo(() => scaleFor("integer_cn", { palette, dark }), [palette, dark]);
  const annotColors = useMemo(() => {
    const m = new Map<string, string>();
    if (!annotation) return m;
    [...new Set(meta.obs[annotation].filter((v) => v != null).map(String))].sort().forEach((k, i) => m.set(k, CATEGORICAL[i % CATEGORICAL.length]));
    return m;
  }, [annotation, meta]);

  const profileData = useMemo(() => {
    if (!cell) return null;
    const L = cell.layers;
    const calls = totalLayer ? L[totalLayer] : null;
    let bins: (number | null)[] | null = null;
    if (binsLayer && L[binsLayer[0]]) {
      const v = L[binsLayer[0]];
      if (binsLayer[1] === "cn" || !calls) bins = v;
      else {
        let s = 0, n = 0;
        for (const x of calls) if (x != null) { s += x; n++; }
        const psi = n ? s / n : 2;
        bins = v.map((x) => (x == null ? null : psi * Math.pow(2, x)));
      }
    }
    const iv = meta.intervals;
    const allele = isAllele && L.major_cn && L.minor_cn
      ? {
          major: L.major_cn,
          minor: L.minor_cn,
          ...(ci && iv.major_cn && iv.minor_cn && L[iv.major_cn[0]]?.some((x) => x != null)
            ? { lo: [L[iv.major_cn[0]], L[iv.minor_cn[0]]] as [typeof L.major_cn, typeof L.major_cn], hi: [L[iv.major_cn[1]], L[iv.minor_cn[1]]] as [typeof L.major_cn, typeof L.major_cn] }
            : {}),
        }
      : null;
    return { calls, bins, allele, reads: countLayer ? L[countLayer] : null };
  }, [cell, totalLayer, binsLayer, countLayer, isAllele, ci, meta]);

  const ready = kind === "heatmap" ? !!heat.matrix : !!profileData;

  const draw = (g: CanvasRenderingContext2D, pxPerPt: number) => {
    const w = widthMm / MM_PER_PT, h = heightMm / MM_PER_PT;
    const common: CommonOpts = { w, h, fontPt, title, subtitle, theme, legend };
    g.setTransform(pxPerPt, 0, 0, pxPerPt, 0, 0);
    if (kind === "heatmap" && heat.matrix) {
      drawHeatmapFigure(g, {
        ...common, meta, matrix: heat.matrix, scale, scaleTitle: isAllele ? "Major | minor copies" : "Copy number",
        shown, rowOf, gx, rows, dendro: dendro && d.dendro ? d.dendro : null, cytobands: ideogram ? cytobands : null,
        annotation: annotation ? { name: annotation.replace(/_/g, " "), values: meta.obs[annotation], colors: annotColors } : null,
        pxPerPt,
      });
    } else if (kind === "profile" && profileData && cell) {
      drawProfileFigure(g, { ...common, meta, cell, gx, ...profileData, reads: reads ? profileData.reads : null, stateScale, colorBins });
    }
  };

  // live preview, fitted to the dialog
  const box = useRef<HTMLDivElement>(null);
  const cv = useRef<HTMLCanvasElement>(null);
  useEffect(() => {
    const c = cv.current, b = box.current;
    if (!c || !b || !ready) return;
    const wPt = widthMm / MM_PER_PT, hPt = heightMm / MM_PER_PT;
    const fitPx = Math.min((b.clientWidth - 24) / wPt, (b.clientHeight - 24) / hPt);
    const dpr = Math.max(1, devicePixelRatio || 1);
    c.width = Math.round(wPt * fitPx * dpr); c.height = Math.round(hPt * fitPx * dpr);
    c.style.width = Math.round(wPt * fitPx) + "px"; c.style.height = Math.round(hPt * fitPx) + "px";
    draw(c.getContext("2d")!, fitPx * dpr);
  });

  const download = async () => {
    setBusy(true);
    await new Promise((r) => setTimeout(r, 30));
    try {
      const W = Math.round((widthMm / 25.4) * dpi), H = Math.round((heightMm / 25.4) * dpi);
      const c = document.createElement("canvas");
      c.width = W; c.height = H;
      draw(c.getContext("2d")!, dpi / 72);
      const base = `${(meta.name || "cnvscope").replace(/[^\w.-]+/g, "_")}_${kind === "heatmap" ? (isAllele ? "allele_heatmap" : "heatmap") : cellName.replace(/[^\w.-]+/g, "_") || "profile"}`;
      if (format === "pdf") {
        const { jsPDF } = await import("jspdf");
        const pdf = new jsPDF({ unit: "mm", format: [widthMm, heightMm], orientation: widthMm > heightMm ? "landscape" : "portrait" });
        pdf.addImage(c.toDataURL("image/png"), "PNG", 0, 0, widthMm, heightMm);
        pdf.save(base + ".pdf");
      } else {
        const blob: Blob = await new Promise((r) => c.toBlob((b) => r(b!), "image/png"));
        const a = document.createElement("a");
        a.href = URL.createObjectURL(blob);
        a.download = `${base}_${dpi}dpi.png`;
        a.click();
        setTimeout(() => URL.revokeObjectURL(a.href), 2000);
      }
    } finally {
      setBusy(false);
    }
  };

  useEffect(() => {
    const h = (e: KeyboardEvent) => { if (e.key === "Escape") { e.stopPropagation(); onClose(); } };
    addEventListener("keydown", h, true);
    return () => removeEventListener("keydown", h, true);
  }, [onClose]);

  const preset = WIDTHS.find(([, v]) => v === widthMm)?.[1] ?? 0;
  const px = `${Math.round((widthMm / 25.4) * dpi)} × ${Math.round((heightMm / 25.4) * dpi)} px`;

  return (
    <div className="modal-back" onMouseDown={(e) => e.target === e.currentTarget && onClose()}>
      <div className="modal" role="dialog" aria-label="Download figure">
        <div className="modal-preview" ref={box}>
          {ready ? <canvas ref={cv} /> : <div className="empty">{kind === "profile" ? "Select a cell first." : "Loading…"}</div>}
        </div>
        <div className="modal-side">
          <div className="modal-head">
            <h2>Download figure</h2>
            <button className="iconbtn" onClick={onClose} aria-label="Close">✕</button>
          </div>
          <div className="modal-opts">
            <Row label="Title"><input className="tinput" value={title} onChange={(e) => setTitle(e.target.value)} /></Row>
            <Row label="Subtitle" hint={subtitle || "none"}><Switch on={showSubtitle} onChange={setShowSubtitle} label="Subtitle" /></Row>
            <Row label="Width">
              <select className="select" value={preset} onChange={(e) => +e.target.value && setWidthMm(+e.target.value)}>
                {WIDTHS.map(([l, v]) => <option key={v} value={v}>{l}</option>)}
                {!preset && <option value={0}>Custom</option>}
              </select>
            </Row>
            <Row label="Size" hint="millimetres">
              <span className="numfield"><input type="number" min={40} max={400} value={widthMm} onChange={(e) => setWidthMm(Math.max(40, Math.min(400, +e.target.value || 40)))} /><span>w</span></span>
              <span className="numfield"><input type="number" min={30} max={400} value={heightMm} onChange={(e) => setHeightMm(Math.max(30, Math.min(400, +e.target.value || 30)))} /><span>h</span></span>
            </Row>
            <Row label="Font size"><Seg value={fontPt} options={[[6, "6"], [7, "7"], [8, "8"], [9, "9"], [10, "10 pt"]]} onChange={setFontPt} /></Row>
            <Row label="Theme"><Seg value={themeName} options={[["print", "Print (light)"], ["screen", "As on screen"]]} onChange={setThemeName} /></Row>
            <Row label="Region" hint={zoomed ? regionLabel(meta, st.zoom!.x) : "zoom the heatmap to export a region"}>
              <Seg value={region} options={zoomed ? [["genome", "Whole genome"], ["zoom", "Zoomed view"]] : [["genome", "Whole genome"]]} onChange={setRegion} />
            </Row>
            <Row label="Legend"><Switch on={legend} onChange={setLegend} label="Legend" /></Row>
            {kind === "heatmap" ? <>
              {d.dendro && <Row label="Dendrogram"><Switch on={dendro} onChange={setDendro} label="Dendrogram" /></Row>}
              {meta.genome && <Row label="Ideogram" hint={cytobands ? meta.genome : "chromosome bands unavailable offline"}><Switch on={ideogram && !!cytobands} onChange={setIdeogram} label="Ideogram" /></Row>}
              {annotCols.length > 0 && (
                <Row label="Annotation track">
                  <select className="select" value={annotation} onChange={(e) => setAnnotation(e.target.value)}>
                    <option value="">none</option>
                    {annotCols.map((c) => <option key={c} value={c}>{c}</option>)}
                  </select>
                </Row>
              )}
            </> : <>
              {countLayer && <Row label="Read depth"><Switch on={reads} onChange={setReads} label="Read depth" /></Row>}
              <Row label="Colour bins by state"><Switch on={colorBins} onChange={setColorBins} label="Colour bins by state" /></Row>
              {isAllele && d.hasCi && <Row label="BAF intervals"><Switch on={ci} onChange={setCi} label="BAF intervals" /></Row>}
            </>}
            <Row label="Format"><Seg value={format} options={[["png", "PNG"], ["pdf", "PDF"]]} onChange={setFormat} /></Row>
            <Row label="Resolution" hint={px}><Seg value={dpi} options={[[300, "300 dpi"], [600, "600 dpi"]]} onChange={setDpi} /></Row>
          </div>
          <div className="modal-foot">
            <span />
            <button className="iconbtn primary" disabled={!ready || busy} onClick={download}>{busy ? "Rendering…" : `Download ${format.toUpperCase()}`}</button>
          </div>
        </div>
      </div>
    </div>
  );
}
