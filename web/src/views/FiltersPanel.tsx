// The Filters panel beside the heatmap: change a QC filter and watch the heatmap react next
// to it. Active filters first; any other metric is one click away. Same filters as Cell QC.

import { useMemo, useState } from "react";
import { cellsTable, defaultFilters, downloadText, isActive } from "../lib/filters";
import { useDerived } from "../state/derived";
import { useStore } from "../state/store";
import { FilterCard, FilterWindow } from "./QCTab";
import { SegmentControl } from "./SegmentControl";
import { DlIcon } from "./widgets";

export function FiltersPanel() {
  const d = useDerived();
  const { meta, vals, metAll, passing, failCounts, ranges, filters, nPass } = d;
  const set = useStore((s) => s.set);
  const setFilter = useStore((s) => s.setFilter);
  const [open, setOpen] = useState<string | null>(null);
  const [adding, setAdding] = useState(false);
  const met = useMemo(() => metAll.filter((m) => ranges[m.key]), [metAll, ranges]);
  const active = met.filter((m) => isActive(filters[m.key]));
  const others = met.filter((m) => !isActive(filters[m.key]));
  // a metric just added stays in the list while its filter is being set
  const [shownExtra, setShownExtra] = useState<string[]>([]);
  const extra = others.filter((m) => shownExtra.includes(m.key));
  const rest = others.filter((m) => !shownExtra.includes(m.key));
  const card = (key: string) => {
    const m = met.find((x) => x.key === key)!;
    return (
      <FilterCard key={key} m={m} values={vals[key]} passing={passing} filter={filters[key]} count={failCounts[key] ?? 0} range={ranges[key]}
        plotted={false} onPlot={() => setOpen(key)} onOpen={() => setOpen(key)} onChange={(f) => setFilter(key, f)} />
    );
  };
  const oi = met.findIndex((m) => m.key === open);
  const go = (k: number) => setOpen(met[(oi + k + met.length) % met.length].key);
  const safe = (meta.name || "cnvscope").replace(/[^\w.-]+/g, "_");

  return (
    <aside className="fpanel" aria-label="QC filters">
      <div className="fpanel-head">
        <div>
          <h2>Filters</h2>
          <div className="ok"><b>{nPass.toLocaleString()}</b> of {meta.n_cells.toLocaleString()} cells pass</div>
        </div>
        <button className="iconbtn" onClick={() => set({ filtersOpen: false })} aria-label="Close the filters panel">✕</button>
      </div>
      <div className="fpanel-bar">
        <button className="iconbtn" onClick={() => set({ filters: defaultFilters(vals, meta.tool) })}>Defaults</button>
        <button className="iconbtn" disabled={!active.length} onClick={() => set({ filters: {} })}>Clear all</button>
        <button className="iconbtn" onClick={() => downloadText(cellsTable(meta, vals, filters, d), `${safe}_cells_qc.tsv`)} title="Every cell with pass/fail, failed filters and all metrics (TSV)"><DlIcon />Cells</button>
      </div>
      <SegmentControl />
      <div className="fpanel-list">
        {active.length === 0 && extra.length === 0 && <p className="ok fpanel-empty">No filters: every cell passes. Add one below.</p>}
        {active.map((m) => card(m.key))}
        {extra.map((m) => card(m.key))}
        <div className="fpanel-add">
          <button className="fbtn" onClick={() => setAdding((a) => !a)} aria-expanded={adding}>{adding ? "Done" : "+ Add a filter"}</button>
          {adding && (
            <div className="fpanel-metrics">
              {rest.map((m) => (
                <button key={m.key} className="fchoice" title={m.help} onClick={() => { setShownExtra((x) => [...x, m.key]); setAdding(false); }}>
                  {m.label}<span>{m.group}</span>
                </button>
              ))}
              {!rest.length && <span className="ok">Every metric is listed.</span>}
            </div>
          )}
        </div>
      </div>
      <p className="hint">Drag a handle or type a value · double-click a histogram to open it large · flagged cells are marked amber in the strip beside the dendrogram</p>
      {oi >= 0 && (
        <FilterWindow m={met[oi]} values={vals[met[oi].key]} passing={passing} filter={filters[met[oi].key]} count={failCounts[met[oi].key] ?? 0}
          range={ranges[met[oi].key]} onChange={(f) => setFilter(met[oi].key, f)} onClose={() => setOpen(null)} onPrev={() => go(-1)} onNext={() => go(1)}
          onPlot={() => { setOpen(null); set({ tab: "qc" }); }} />
      )}
    </aside>
  );
}
