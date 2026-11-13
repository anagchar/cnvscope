import { useEffect, useRef } from "react";
import type { Meta } from "./api/client";
import { defaultFilters } from "./lib/filters";
import { DerivedProvider, useDerived } from "./state/derived";
import { useStore } from "./state/store";
import { DatasetTab } from "./views/DatasetTab";
import { FiltersPanel } from "./views/FiltersPanel";
import { HeatmapCard } from "./views/HeatmapCard";
import { CellQCCard, LorenzCard, ProfileCard } from "./views/ProfileCard";
import { QCTab } from "./views/QCTab";
import { Sidebar, navigate } from "./views/Sidebar";
import { LogoIcon } from "./views/Logo";
import { UploadScreen } from "./views/UploadScreen";
import { TAB_ICONS, ThemeIcon } from "./views/widgets";

export function App() {
  const meta = useStore((s) => s.meta);
  const booting = useStore((s) => s.booting);
  const picking = useStore((s) => s.picking);
  const boot = useStore((s) => s.boot);
  useEffect(() => { boot(); }, [boot]);
  if (booting) return null;
  if (!meta || picking) return <UploadScreen />;
  return (
    <DerivedProvider meta={meta} key={meta.name + meta.n_cells}>
      <Explorer meta={meta} />
    </DerivedProvider>
  );
}

function Explorer({ meta }: { meta: Meta }) {
  const d = useDerived();
  const s = useStore();
  const hasCounts = !!d.countLayer;
  const searchRef = useRef<HTMLInputElement>(null);

  // defaults for a new dataset: QC filters from the cells' own distribution
  useEffect(() => {
    const st = useStore.getState();
    if (st.filters == null) st.set({ filters: defaultFilters(d.vals, meta.tool) });
    // the heatmap's track: the tool's clones, else a cell-type column (e.g. from --annotations)
    const track = ["clone", "cell_type"].find((k) => meta.obs_kinds[k] === "categorical");
    if (st.annotation == null && track) st.set({ annotation: track });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [meta]);
  useEffect(() => {
    const st = useStore.getState();
    if (st.selected == null && d.shown.length) st.set({ selected: d.shown[0] });
  }, [d.shown]);
  // zoom resets when the set of cells shown changes
  const nShown = d.shown.length;
  useEffect(() => { useStore.getState().set({ zoom: null }); }, [nShown]);

  useEffect(() => {
    const h = (e: KeyboardEvent) => {
      const t = e.target as HTMLElement;
      if (t.tagName === "INPUT" || t.tagName === "SELECT") { if (e.key === "Escape") t.blur(); return; }
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      if (e.key === "ArrowUp" || e.key === "k") { e.preventDefault(); navigate(-1); }
      else if (e.key === "ArrowDown" || e.key === "j") { e.preventDefault(); navigate(1); }
      else if (e.key === "Escape") useStore.getState().set({ zoom: null });
      else if (e.key === "/") { e.preventDefault(); useStore.getState().set({ panel: true }); setTimeout(() => searchRef.current?.focus(), 0); }
      else if (e.key === "t") useStore.getState().toggleTheme();
    };
    addEventListener("keydown", h);
    return () => removeEventListener("keydown", h);
  }, []);

  const profile = <ProfileCard />;
  return (
    <div className={"app" + (s.panel ? "" : " nopanel")}>
      <header>
        <div className="brand">
          <LogoIcon size={36} />
          cnvscope
        </div>
        <div className="project" title={meta.name}>{meta.name}</div>
        <div className="chips">
          <span className="chip">{meta.n_cells.toLocaleString()} cells · {meta.n_bins.toLocaleString()} bins</span>
          <span className="chip">{meta.modality} · {meta.tool}</span>
          {d.hasAS && <span className="chip as">allele-specific</span>}
        </div>
        <div className="spacer" />
        <button className="iconbtn" onClick={() => s.set({ panel: !s.panel })} title={s.panel ? "Hide the side panel" : "Show the side panel"} aria-pressed={s.panel}>
          <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5"><rect x="1.5" y="2.5" width="13" height="11" rx="2" /><path d="M6 2.5v11" /></svg><span className="btnlabel">Panel</span>
        </button>
        <button className="iconbtn" onClick={() => s.set({ picking: true, error: null })} title="Open another file">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5-5 5 5M12 5v12" /></svg><span className="btnlabel">Open</span>
        </button>
        <button className="iconbtn" onClick={s.toggleTheme} title="Light or dark (t)" aria-label="Switch between light and dark"><ThemeIcon /></button>
      </header>

      <Sidebar ref={searchRef} />

      <main>
        <div className="tabs" role="tablist">
          {([["cn", "Copy number", "Heatmap and cell profiles"], ["qc", "Cell QC", "Per-cell quality metrics"], ["data", "Dataset", "What this dataset holds"]] as const).map(([id, label, title]) => (
            <button key={id} className={"tab" + (s.tab === id ? " on" : "")} role="tab" aria-selected={s.tab === id} title={title} onClick={() => s.set({ tab: id })}>{TAB_ICONS[id]}{label}</button>
          ))}
        </div>
        {s.error && <div className="error">{s.error}</div>}
        {s.tab === "cn" && <section className={"cn-layout" + (s.filtersOpen ? " with-filters" : "")}>
          <div className="cn-main">

          <HeatmapCard />
          {profile}
          {s.selected != null && (hasCounts ? <div className="grid2"><LorenzCard /><CellQCCard /></div> : <CellQCCard />)}
          </div>
          {s.filtersOpen && <FiltersPanel />}
        </section>}
        {s.tab === "qc" && <section><QCTab>{profile}</QCTab></section>}
        {s.tab === "data" && <section><DatasetTab /></section>}
      </main>
    </div>
  );
}
