import { defineConfig } from "vitest/config";
import react from "@vitejs/plugin-react";

// The built UI ships inside the Python wheel; `cnvscope serve` serves it.
// `vite build --mode pages` (npm run build:pages) makes the GitHub Pages site
// instead: no server, files are read in the browser (.env.pages: VITE_STATIC=1).
export default defineConfig(({ mode }) => ({
  plugins: [react()],
  base: "./",
  build: { outDir: mode === "pages" ? "dist-pages" : "../src/cnvscope/static", emptyOutDir: true },
  server: { proxy: { "/api": "http://127.0.0.1:8765" } },
  test: { environment: "node" },
}));
