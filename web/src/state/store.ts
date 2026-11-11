import { create } from "zustand";
import { fetchLayer, fetchMeta, openDemo, openFile, type Meta } from "../api/client";
import type { Filter, Filters } from "../lib/filters";
import type { PaletteName } from "../lib/palettes";

export type Theme = "light" | "dark";
export type Tab = "cn" | "qc" | "data";
/** Genomic window (cumulative coordinates) and display-row window of the heatmap. */
export interface Zoom {
  x: [number, number];
  y: [number, number];
}
/** The heatmap's track: a layer name, or ALLELE for major|minor states. */
export const ALLELE = "__allele__";

interface State {
  meta: Meta | null;
  booting: boolean;
  uploading: boolean;
  /** show the upload screen even though a dataset is open ("Open" in the header) */
  picking: boolean;
  error: string | null;
  /** loaded layers, rows in dataset order */
  matrices: Record<string, Float32Array>;
  loading: Record<string, number>;

  view: string;
  selected: number | null;
  zoom: Zoom | null;
  tab: Tab;
  panel: boolean;
  theme: Theme;
  palette: PaletteName;
  showDendro: boolean;
  showCi: boolean;
  passOnly: boolean;
  minSegmentMb: number;
  heatmapH: number;
  /** QC filters by metric; null until the dataset's defaults are set */
  filters: Filters | null;
  /** the Filters panel beside the heatmap */
  filtersOpen: boolean;
  cellType: string;
  annotation: string | null;
  search: string;

  boot: () => Promise<void>;
  open: (file: File) => Promise<void>;
  demo: () => Promise<void>;
  setMeta: (m: Meta) => void;
  ensureLayer: (name: string) => Promise<void>;
  setView: (v: string) => void;
  set: (p: Partial<State>) => void;
  toggleTheme: () => void;
  setPalette: (p: PaletteName) => void;
  setFilter: (metric: string, f: Filter | null) => void;
}

// one palette for everyone: the SIGNALS / scgenome colours (lib/palettes.ts)
const initialPalette = (): PaletteName => "signals";
const initialTheme = (): Theme => (document.documentElement.dataset.theme === "dark" ? "dark" : "light");

/** Layers the heatmap needs for a view. */
export function layersFor(view: string): string[] {
  return view === ALLELE ? ["major_cn", "minor_cn"] : [view];
}

export const useStore = create<State>((set, get) => ({
  meta: null,
  booting: true,
  uploading: false,
  picking: false,
  error: null,
  matrices: {},
  loading: {},
  view: "",
  selected: null,
  zoom: null,
  tab: "cn",
  panel: true,
  theme: initialTheme(),
  palette: initialPalette(),
  showDendro: true,
  showCi: true,
  passOnly: true,
  minSegmentMb: 0,
  heatmapH: 560,
  filters: null,
  filtersOpen: false,
  cellType: "All",
  annotation: null,
  search: "",

  boot: async () => {
    try {
      const m = await fetchMeta();
      if (m) get().setMeta(m);
    } catch (e) {
      set({ error: String(e instanceof Error ? e.message : e) });
    } finally {
      set({ booting: false });
    }
  },

  open: async (file) => {
    set({ uploading: true, error: null });
    try {
      get().setMeta(await openFile(file));
    } catch (e) {
      set({ error: e instanceof Error ? e.message : String(e) });
    } finally {
      set({ uploading: false });
    }
  },

  demo: async () => {
    set({ uploading: true, error: null });
    try {
      get().setMeta(await openDemo());
    } catch (e) {
      set({ error: e instanceof Error ? e.message : String(e) });
    } finally {
      set({ uploading: false });
    }
  },

  setMeta: (meta) => {
    // a new dataset: every view setting goes back to its default
    set({
      meta,
      picking: false,
      error: null,
      matrices: {},
      loading: {},
      view: meta.default_layer,
      selected: null,
      zoom: null,
      tab: "cn",
      showDendro: true,
      minSegmentMb: 0,
      passOnly: true,
      cellType: "All",
      search: "",
      filters: null,
      annotation: null,
    });
    get().ensureLayer(meta.default_layer);
  },

  ensureLayer: async (name) => {
    const { meta, matrices, loading } = get();
    if (!meta || matrices[name] || name in loading) return;
    set((s) => ({ loading: { ...s.loading, [name]: 0 } }));
    try {
      const m = await fetchLayer(name, meta.n_cells, meta.n_bins, (rows) =>
        set((s) => ({ loading: { ...s.loading, [name]: rows } })),
      );
      if (get().meta !== meta) return;
      set((s) => {
        const loading = { ...s.loading };
        delete loading[name];
        return { matrices: { ...s.matrices, [name]: m }, loading };
      });
    } catch (e) {
      set({ error: e instanceof Error ? e.message : String(e) });
    }
  },

  setView: (view) => {
    const meta = get().meta;
    if (!meta) return;
    set({ view });
    for (const l of layersFor(view)) get().ensureLayer(l);
  },

  set: (p) => set(p),

  toggleTheme: () => {
    const theme: Theme = get().theme === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = theme;
    try {
      localStorage.setItem("cnvscope-theme", theme);
    } catch {
      /* storage blocked */
    }
    set({ theme });
  },

  setFilter: (metric, f) => {
    const filters = { ...(get().filters ?? {}) };
    if (f && (f.min != null || f.max != null)) filters[metric] = f;
    else delete filters[metric];
    set({ filters });
  },

  setPalette: (palette) => {
    try {
      localStorage.setItem("cnvscope-palette", palette);
    } catch {
      /* storage blocked */
    }
    set({ palette });
  },
}));
