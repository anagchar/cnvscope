import { forwardRef, useEffect, useMemo } from "react";
import { fmt } from "../lib/canvas";
import { cellsTable, describe, downloadText, isActive } from "../lib/filters";
import { infoFor } from "../lib/metrics";
import { useDerived, type Derived } from "../state/derived";
import { ALLELE, useStore } from "../state/store";
import { DlIcon, Slider, Switch } from "./widgets";

const LIST_CAP = 400;
const SEM_LABEL: Record<string, string> = {
  integer_cn: "copy number", allele_cn: "allele CN", ternary_state: "gain / loss", continuous_ratio: "log ratio",
  continuous_cn: "CN scale", ratio: "ratio", baf: "BAF", probability: "probability", count: "counts",
};

/** Open Cell QC at the filters (optionally at one metric's card). */
export function editFilters(metric?: string) {
  const st = useStore.getState();
  if (st.tab === "cn") {
    // beside the heatmap, so its reaction shows at once
    st.set({ filtersOpen: true });
    setTimeout(() => {
      const el = metric ? document.querySelector(`.fpanel #filter-${metric}`) : null;
      el?.scrollIntoView({ behavior: "smooth", block: "center" });
      if (el) { el.classList.add("flash"); setTimeout(() => el.classList.remove("flash"), 1200); }
    }, 60);
    return;
  }
  st.set({ tab: "qc" });
  setTimeout(() => {
    const el = document.getElementById(metric ? `filter-${metric}` : "qc-filters");
    el?.scrollIntoView({ behavior: "smooth", block: metric ? "center" : "start" });
    if (metric && el) { el.classList.add("flash"); setTimeout(() => el.classList.remove("flash"), 1200); }
  }, 60);
}

/** Download every cell with pass/fail, the failed filters and all metrics. */
export function downloadCells(d: Derived) {
  const s = useStore.getState();
  const ev = { passing: d.passing, reasons: d.reasons, failCounts: d.failCounts, nPass: d.nPass };
  const safe = (d.meta.name || "cnvscope").replace(/[^\w.-]+/g, "_");
  downloadText(cellsTable(d.meta, d.vals, s.filters ?? {}, ev), `${safe}_cells_qc.tsv`);
}

export function navigate(dir: number) {
  const { selected, set } = useStore.getState();
  const shown = currentShown;
  if (!shown.length) return;
  const i = selected == null ? -1 : shown.indexOf(selected);
  set({ selected: shown[Math.max(0, Math.min(shown.length - 1, i < 0 ? 0 : i + dir))] });
}
// the cells shown, for keyboard navigation outside React
let currentShown: number[] = [];

export const Sidebar = forwardRef<HTMLInputElement>(function Sidebar(_, searchRef) {
  const d = useDerived();
  const { meta, hasAS, hasCi, cellTypes, failCounts, passing, reasons, nPass, shown, vals, sem } = d;
  const activeFilters = Object.entries(d.filters).filter(([, f]) => isActive(f));
  const s = useStore();
  currentShown = shown;

  const isAllele = s.view === ALLELE;
  const totalView = meta.default_layer;
  // scDNA shows Total and, when there is one, Allele-specific; the per-bin tracks (RDR,
  // read counts, ...) stay in the cell profile. Other modalities list every layer.
  const others = meta.modality === "scDNA"
    ? []
    : meta.layers.filter((l) => l.name !== totalView && !(hasAS && (l.name === "major_cn" || l.name === "minor_cn")));
  const annotCols = Object.keys(meta.obs_kinds).filter((c) => meta.obs_kinds[c] === "categorical");

  const list = useMemo(() => {
    const q = s.search.trim().toLowerCase();
    return q ? shown.filter((i) => meta.cells[i].toLowerCase().includes(q)) : shown;
  }, [shown, s.search, meta]);
  useEffect(() => { document.querySelector(".cell.on")?.scrollIntoView({ block: "nearest" }); }, [s.selected]);

  return (
    <aside><div className="aside-scroll">
      <div className="group">
        <h3>Copy number</h3>
        {hasAS && (
          <div className="seg">
            <button className={isAllele ? "" : s.view === totalView ? "on" : ""} onClick={() => s.setView(totalView)}>Total</button>
            <button className={isAllele ? "on" : ""} onClick={() => s.setView(ALLELE)}>Allele-specific</button>
          </div>
        )}
        {others.length > 0 && (
          <div className="layers">
            {!hasAS && <button className={s.view === totalView ? "on" : ""} onClick={() => s.setView(totalView)}><span>{totalView}</span><small>{SEM_LABEL[sem[totalView]]}</small></button>}
            {others.map((l) => (
              <button key={l.name} className={s.view === l.name ? "on" : ""} onClick={() => s.setView(l.name)}>
                <span>{l.name}</span><small>{SEM_LABEL[l.semantics] ?? l.semantics}</small>
              </button>
            ))}
          </div>
        )}
        {d.dendro && <div className="row"><div><label>Dendrogram</label><small>Ward clustering of profiles</small></div><Switch on={s.showDendro} onChange={(v) => s.set({ showDendro: v })} label="Dendrogram" /></div>}
        {isAllele && hasCi && <div className="row"><div><label>BAF intervals</label><small>5–95% ribbons on the profile</small></div><Switch on={s.showCi} onChange={(v) => s.set({ showCi: v })} label="BAF intervals" /></div>}
        <Slider label="Heatmap height" value={s.heatmapH} shown={`${s.heatmapH} px`} min={250} max={900} step={10} onChange={(v) => s.set({ heatmapH: v })} />
        {annotCols.length > 0 && (
          <div className="row">
            <div><label>Annotation track</label><small>a colour strip beside the heatmap</small></div>
            <select className="select" value={s.annotation ?? ""} onChange={(e) => s.set({ annotation: e.target.value || null })}>
              <option value="">none</option>
              {annotCols.map((c) => <option key={c} value={c}>{c}</option>)}
            </select>
          </div>
        )}
      </div>

      <div className="group">
        <h3>Cells</h3>
        <div className="qcsum">
          <div className="qcsum-n"><b>{nPass.toLocaleString()}</b> of {meta.n_cells.toLocaleString()} cells pass QC</div>
          {activeFilters.length > 0 ? (
            <div className="chips-wrap">
              {activeFilters.map(([k, f]) => (
                <span key={k} className="fchip" title={`${failCounts[k] ?? 0} cells flagged`}>
                  <button className="fchip-l" onClick={() => editFilters(k)}>{infoFor(k).label} {describe(k, f)}</button>
                  <button className="fchip-x" aria-label={`Remove the ${infoFor(k).label} filter`} onClick={() => s.setFilter(k, null)}>✕</button>
                </span>
              ))}
            </div>
          ) : <small className="ok">No QC filters: every cell passes.</small>}
          <div className="qcsum-btns">
            <button className="iconbtn" onClick={() => editFilters()}>Edit filters</button>
            <button className="iconbtn" onClick={() => downloadCells(d)} title="Every cell with pass/fail, failed filters and all metrics (TSV)"><DlIcon />Cells</button>
          </div>
        </div>
        <div className="row">
          <div><label>Passing QC only</label><small>{nPass === meta.n_cells ? "every cell passes the filters" : `hides ${(meta.n_cells - nPass).toLocaleString()} flagged cells`}</small></div>
          <Switch on={s.passOnly} onChange={(v) => s.set({ passOnly: v })} label="Passing QC only" />
        </div>
        {cellTypes.length > 0 && (
          <div className="row" style={{ display: "block" }}>
            <label>Cell type</label>
            <div className="seg wrap" style={{ marginTop: 8 }}>
              {["All", ...cellTypes].map((ct) => <button key={ct} className={s.cellType === ct ? "on" : ""} onClick={() => s.set({ cellType: ct })}>{ct}</button>)}
            </div>
          </div>
        )}
        <div className="row" style={{ background: "transparent", padding: "6px 2px 0" }}><span className="ok">{shown.length.toLocaleString()} cells shown</span></div>
      </div>

      <div className="group">
        <h3>Find a cell</h3>
        <input ref={searchRef} className="search" value={s.search} onChange={(e) => s.set({ search: e.target.value })} placeholder="Barcode or name  ( / )" autoComplete="off" spellCheck={false} />
        <div className="cells" role="listbox" aria-label="Cells">
          {list.slice(0, LIST_CAP).map((i) => {
            const ok = passing[i];
            const ct = meta.obs.cell_type?.[i];
            return (
              <button key={i} className={"cell" + (i === s.selected ? " on" : "")} role="option" aria-selected={i === s.selected} onClick={() => s.set({ selected: i })}>
                <span className="dot" style={{ background: ok ? "var(--good)" : "var(--serious)" }} title={ok ? "passes QC" : "flagged: " + reasons[i]} />
                <span className="nm">{meta.cells[i]}</span>
                {ct != null && <span className="ct">{String(ct)}</span>}
                <span className="pl">{vals.ploidy ? fmt(vals.ploidy[i]) : ""}</span>
              </button>
            );
          })}
          {list.length > LIST_CAP && <div className="more">and {(list.length - LIST_CAP).toLocaleString()} more; refine the search</div>}
          {!list.length && <div className="more">No match.</div>}
        </div>
        <div className="navbtns"><button className="iconbtn" onClick={() => navigate(-1)}>↑ Previous</button><button className="iconbtn" onClick={() => navigate(1)}>Next ↓</button></div>
      </div>
      <div className="hint-keys"><kbd>↑</kbd> <kbd>↓</kbd> cells · <kbd>drag</kbd> zoom · <kbd>Esc</kbd> reset · <kbd>/</kbd> search · <kbd>t</kbd> theme</div>
    </div></aside>
  );
});
