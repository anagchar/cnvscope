import { useRef, useState } from "react";
import { isBrowserOnly } from "../api/client";
import { useStore } from "../state/store";
import { LogoIcon } from "./Logo";
import { ThemeIcon } from "./widgets";

export function UploadScreen() {
  const { open, demo, uploading, error, meta, toggleTheme, set } = useStore();
  const [over, setOver] = useState(false);
  const file = useRef<HTMLInputElement>(null);
  const take = (f?: File) => { if (f) open(f); };
  const local = isBrowserOnly();
  return (
    <div className="upload">
      <button className="iconbtn" style={{ position: "fixed", top: 16, right: 16 }} onClick={toggleTheme} title="Light or dark (t)" aria-label="Switch between light and dark"><ThemeIcon /></button>
      <div className="inner">
        <div className="hero">
          <div className="mark"><LogoIcon size={150} /></div>
          <h1>cnvscope</h1>
          <p>Single-cell copy number, interactively</p>
        </div>
        <button
          className={"card drop" + (over ? " over" : "")}
          onDragOver={(e) => { e.preventDefault(); setOver(true); }}
          onDragLeave={() => setOver(false)}
          onDrop={(e) => { e.preventDefault(); setOver(false); take(e.dataTransfer.files[0]); }}
          onClick={() => file.current?.click()}
        >
          <input ref={file} type="file" accept={local ? ".json,.gz,.h5ad" : ".json,.gz,.h5ad,.rds"} hidden onChange={(e) => take(e.target.files?.[0])} />
          {uploading ? <><div className="spinner" /><div className="small">Opening…</div></> : <>
            <svg width="34" height="34" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5-5 5 5M12 5v12" /></svg>
            <div className="big">Drop <code>ascat_data.json</code> or a <code>.h5ad</code> here</div>
            <div className="small">or click to choose a file · .json, .json.gz or .h5ad{local ? " · opened in your browser, nothing is uploaded" : ""}</div>
          </>}
        </button>
        {error && <div className="err">{error}</div>}
        <div className="actions">
          <button className="iconbtn" onClick={demo} disabled={uploading}>Load the demo dataset · 300 cells, allele-specific</button>
          {meta && <button className="iconbtn" onClick={() => set({ picking: false, error: null })}>Back to {meta.name}</button>}
        </div>
        <div className="card">
          <div className="hd"><h2>From the command line</h2><span className="sub">{local ? "read ASCAT.sc, HiScanner, CHISEL or epiAneufinder output into a .h5ad, then open it here" : "read a tool's output, save it as .h5ad, list the readers"}</span></div>
          <div className="bd"><pre className="mono">{`${local ? "pip install git+https://github.com/anagchar/cnvscope\n" : ""}cnvscope --mode dna --tool chisel --input chisel_out/ --outdir results/
cnvscope --mode atac --epianeufinder --input epi_out/ --outdir results/
cnvscope tools`}</pre></div>
        </div>
        <div className="card">
          <div className="hd"><h2>From an ASCAT.sc pipeline (R)</h2><span className="sub">writes a small viewer file beside the result</span></div>
          <div className="bd"><pre className="mono">{`source(system("cnvscope rscript", intern = TRUE))\nexport_cnvscope(res, "results/result_object.Rda")`}</pre></div>
        </div>
      </div>
    </div>
  );
}
