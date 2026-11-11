import { useEffect, useRef, useState, type CSSProperties, type RefObject } from "react";
import type { Threshold } from "../lib/metrics";

/** Width of an element, followed as it resizes. */
export function useWidth(ref: RefObject<HTMLElement | null>, fallback = 600): number {
  const [w, setW] = useState(fallback);
  useEffect(() => {
    if (!ref.current) return;
    const ro = new ResizeObserver((e) => {
      const v = e[0]?.contentRect.width;
      if (v > 0) setW(Math.round(v));
    });
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, [ref]);
  return w;
}

/** Width and height of a `.fill` box: charts that take their card's full height. */
export function useSize(ref: RefObject<HTMLElement | null>, minH: number) {
  const [s, setS] = useState({ w: 600, h: minH });
  useEffect(() => {
    if (!ref.current) return;
    const ro = new ResizeObserver((e) => {
      const r = e[0]?.contentRect;
      if (!r || !(r.width > 0)) return;
      setS({ w: Math.round(r.width), h: Math.max(minH, Math.round(r.height)) });
    });
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, [ref, minH]);
  return s;
}

export function ThemeIcon() {
  return (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5">
      <circle cx="8" cy="8" r="6" />
      <path d="M8 2a6 6 0 0 0 0 12z" fill="currentColor" />
    </svg>
  );
}

export function DlIcon() {
  return (
    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5 5 5-5M12 15V3" />
    </svg>
  );
}

export function Switch({ on, onChange, label }: { on: boolean; onChange: (v: boolean) => void; label: string }) {
  return <button className="switch" role="switch" aria-checked={on} aria-label={label} onClick={() => onChange(!on)} />;
}

export function DownloadMenu({ onSave }: { onSave: (f: "png" | "pdf") => void }) {
  const [open, setOpen] = useState(false);
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!open) return;
    const h = (e: MouseEvent) => {
      if (!ref.current?.contains(e.target as Node)) setOpen(false);
    };
    document.addEventListener("mousedown", h);
    return () => document.removeEventListener("mousedown", h);
  }, [open]);
  return (
    <div className="menu" ref={ref}>
      <button className="iconbtn" onClick={() => setOpen((o) => !o)} title="Save the figure" aria-haspopup="menu" aria-expanded={open}>
        <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
          <path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5 5 5-5M12 15V3" />
        </svg>
        Save
      </button>
      {open && (
        <div className="pop" role="menu">
          {(["png", "pdf"] as const).map((f) => (
            <button key={f} role="menuitem" onClick={() => { onSave(f); setOpen(false); }}>
              {f.toUpperCase()}
            </button>
          ))}
        </div>
      )}
    </div>
  );
}

/** A slider whose far end means "off": cells beyond the value are flagged. */
export function ThresholdSlider({ t, value, range, count, passOnly, onChange }: {
  t: Threshold;
  value: number | null;
  range: [number, number];
  count: number;
  passOnly: boolean;
  onChange: (v: number | null) => void;
}) {
  const [lo, hi] = range, span = hi - lo || 1;
  const off = value == null;
  const pos = off ? (t.dir === "max" ? 1 : 0) : Math.max(0, Math.min(1, (value - lo) / span));
  return (
    <div className="slider">
      <div className="top"><span>{t.label}</span><span>{off ? "off" : t.f(value)}</span></div>
      <input
        type="range" min="0" max="1" step="0.001" value={pos} aria-label={t.label}
        style={{ "--p": (pos * 100).toFixed(1) + "%" } as CSSProperties}
        onChange={(e) => {
          const p = +e.target.value, isOff = t.dir === "max" ? p >= 0.999 : p <= 0.001, v = lo + p * span;
          onChange(isOff ? null : t.int ? Math.round(v) : v);
        }}
      />
      <small className={"hint" + (count ? " on" : "")}>
        {count ? `${count.toLocaleString()} cell${count > 1 ? "s" : ""} ${passOnly ? "hidden" : "flagged"} by this threshold` : t.hint}
      </small>
    </div>
  );
}

export function Slider({ label, value, shown, min, max, step, hint, onChange }: {
  label: string; value: number; shown: string; min: number; max: number; step: number; hint?: string; onChange: (v: number) => void;
}) {
  const p = ((100 * (value - min)) / (max - min || 1)).toFixed(1) + "%";
  return (
    <div className="slider">
      <div className="top"><span>{label}</span><span>{shown}</span></div>
      <input type="range" min={min} max={max} step={step} value={value} aria-label={label} style={{ "--p": p } as CSSProperties} onChange={(e) => onChange(+e.target.value)} />
      {hint && <small className="hint">{hint}</small>}
    </div>
  );
}

export const TAB_ICONS = {
  cn: (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
      <path d="M1.5 9.5h3v-4h3v6h3v-8h2.5" />
      <path d="M1.5 14h13" opacity=".45" />
    </svg>
  ),
  qc: (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
      <circle cx="4" cy="11" r="1.3" /><circle cx="7.5" cy="6.5" r="1.3" /><circle cx="11" cy="9" r="1.3" /><circle cx="13" cy="4" r="1.3" />
      <path d="M1.5 14.5h13M1.5 1.5v13" opacity=".45" />
    </svg>
  ),
  data: (
    <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
      <rect x="2.5" y="1.5" width="11" height="13" rx="2" />
      <path d="M5.5 5.5h5M5.5 8.5h5M5.5 11.5h3" />
    </svg>
  ),
};
