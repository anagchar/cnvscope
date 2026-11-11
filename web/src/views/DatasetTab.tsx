import { fmtBp, fmtCount, quantile } from "../lib/canvas";
import { useDerived } from "../state/derived";

export function DatasetTab() {
  const { meta, vals, hasAS, hasCi, cellTypes, binSizeMb, metAll, dendro } = useDerived();
  const nb = meta.n_bins;
  const reads = vals.n_reads_total ? Array.from(vals.n_reads_total).filter((v) => !Number.isNaN(v)) : null;
  const tiles: [string, string | number, string][] = [
    ["Modality", meta.modality, meta.tool],
    ["Cells", meta.n_cells.toLocaleString(), cellTypes.length ? `${cellTypes.length} cell types` : ""],
  ];
  if (reads?.length) tiles.push(["Reads", fmtCount(reads.reduce((a, v) => a + v, 0)), `median ${fmtCount(quantile([...reads].sort((a, b) => a - b), 0.5))} per cell`]);
  tiles.push(["Bins", nb.toLocaleString(), binSizeMb ? fmtBp(binSizeMb * 1e6) + " each" : ""]);
  tiles.push(["Chromosomes", meta.chr_info.length, fmtBp(meta.bins.end_cum[nb - 1]) + " in bins" + (meta.genome ? ` · ${meta.genome}` : "")]);
  tiles.push(["Allele-specific", hasAS ? "yes" : "no", hasAS ? (hasCi ? "with BAF intervals" : "major and minor copies") : "total copy number only"]);
  tiles.push(["Clustering", dendro ? "Ward.D2" : "–", dendro ? "Manhattan distance of profiles" : "no dendrogram"]);
  return (
    <>
      <div className="card">
        <div className="hd"><h2>Dataset</h2><span className="sub">{meta.name}</span></div>
        <div className="bd" style={{ paddingBottom: 0 }}><div className="tiles">
          {tiles.map(([k, v, s]) => (
            <div className="tile" key={k}><div className="k">{k}</div><div className={"v" + (String(v).length > 9 ? " long" : "")} title={String(v)}>{v}</div><div className="s">{s}</div></div>
          ))}
        </div></div>
      </div>
      <div className="card">
        <div className="hd"><h2>Layers</h2><span className="sub">the per-bin tracks of every cell, and what their values mean</span></div>
        <div className="bd tablewrap"><table>
          <thead><tr><th>layer</th><th>meaning</th><th>interval</th></tr></thead>
          <tbody>{meta.layers.map((l) => (
            <tr key={l.name}><td className="mono">{l.name}</td><td>{l.semantics.replace(/_/g, " ")}</td><td className="ok">{meta.intervals[l.name] ? meta.intervals[l.name].join(" – ") : ""}</td></tr>
          ))}</tbody>
        </table></div>
      </div>
      <div className="card">
        <div className="hd"><h2>QC metrics</h2><span className="sub">from the tool's own export where given, the rest computed as ASCAT.sc's getQCs() does</span></div>
        <div className="bd"><span className="ok">{metAll.map((m) => m.label).join(" · ") || "none"}</span></div>
      </div>
    </>
  );
}
