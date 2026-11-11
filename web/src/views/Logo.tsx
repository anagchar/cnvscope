// The cnvscope mark: the ASCAT.sc karyotype (four chromosomes, centromere dots, red gains,
// a blue loss), drawn as in the original: chromosomes filled with the background colour and
// outlined in the text colour, so they are white with a dark outline on a light page and dark
// with a white outline on a dark page. No tile. Colours: --logo-* in theme.css.
// Downloadable copies: web/public/brand, made by scripts/make_logos.py from the same geometry.

const W = 9, SW = 1.8, CEN = 57, INSET = 1.9;
/** [x, top, height, bands: [top, height, kind]] */
export const CHROMOSOMES: [number, number, number, [number, number, "gain" | "loss"][]][] = [
  [22.5, 28, 64, [[37.6, 12.8, "gain"], [63.2, 9.6, "gain"]]],
  [44.5, 35, 50, [[45, 17.5, "loss"]]],
  [66.5, 25, 70, [[53, 14, "gain"]]],
  [88.5, 38, 44, []],
];
/** the mark's own box, outline included: x 22.5–97.5, y 25–95 */
const VIEW = "21 24 78 72";

export function LogoIcon({ size = 32, title = "cnvscope" }: { size?: number; title?: string }) {
  return (
    <svg className="logo-icon" height={size} width={(size * 78) / 72} viewBox={VIEW} role="img" aria-label={title}>
      {CHROMOSOMES.map(([x, y, h, bands], i) => (
        <g key={i}>
          <rect x={x} y={y} width={W} height={h} rx={W / 2} fill="var(--logo-fill)" stroke="var(--logo-line)" strokeWidth={SW} />
          {bands.map(([by, bh, k], j) => (
            <rect key={j} x={x + INSET} y={by} width={W - 2 * INSET} height={bh} rx={1.6} fill={`var(--logo-${k})`} />
          ))}
          <circle cx={x + W / 2} cy={CEN} r={1.7} fill="var(--logo-line)" />
        </g>
      ))}
    </svg>
  );
}
