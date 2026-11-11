import type { CSSProperties } from "react";
import { isCall } from "../lib/scales";
import { useDerived } from "../state/derived";
import { ALLELE, useStore } from "../state/store";

/**
 * Minimum segment size: runs of a non-baseline call shorter than this show as the baseline
 * (copy number 2) in the heatmap, and "Bins with CNA" (and its filter) follow.
 */
export function SegmentControl({ compact = false }: { compact?: boolean }) {
  const { meta, binSizeMb, sem, totalLayer } = useDerived();
  const view = useStore((s) => s.view);
  const mb = useStore((s) => s.minSegmentMb);
  const set = useStore((s) => s.set);
  const viewSem = view === ALLELE ? "allele_cn" : sem[view];
  if (!totalLayer && !(viewSem && isCall(viewSem))) return null;
  // up to 200 Mb (whole chromosome arms): finer control where it matters
  const max = Math.max(1, Math.min(200, Math.round((meta.n_bins * binSizeMb) / 2)));
  const step = Math.max(1, Math.ceil(binSizeMb));
  const p = ((100 * mb) / max).toFixed(1) + "%";
  return (
    <div className={"segctl" + (compact ? " compact" : "") + (mb ? " on" : "")}>
      <div className="segctl-top">
        <label htmlFor="segctl">Min segment size</label>
        <span>{mb ? `${mb} Mb` : "off"}</span>
      </div>
      <input id="segctl" type="range" min={0} max={max} step={step} value={mb} aria-label="Minimum segment size"
        style={{ "--p": p } as CSSProperties} onChange={(e) => set({ minSegmentMb: +e.target.value })} />
      <small>{viewSem === "ternary_state" ? "Shorter gains and losses show as normal." : "Shorter gains and losses show as copy number 2; Bins with CNA follows."}</small>
    </div>
  );
}
