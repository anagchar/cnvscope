# ruff: noqa: E501 (SVG markup)
"""Write the cnvscope logo files (web/public/brand and the favicon) from one geometry.

The ASCAT.sc karyotype (four chromosomes, centromere dots, red gains, a blue loss) as in the
original: chromosomes filled with the background colour and outlined, no tile. Same geometry as
web/src/views/Logo.tsx. Run: python scripts/make_logos.py
"""

from pathlib import Path

OUT = Path(__file__).parents[1] / "web" / "public" / "brand"
W, SW, CEN, INSET = 9, 1.8, 57, 1.9
CHROMOSOMES = [
    (22.5, 28, 64, [(37.6, 12.8, "gain"), (63.2, 9.6, "gain")]),
    (44.5, 35, 50, [(45, 17.5, "loss")]),
    (66.5, 25, 70, [(53, 14, "gain")]),
    (88.5, 38, 44, []),
]
VIEW = (21, 24, 78, 72)
INK, PAPER, RED, BLUE = "#231F20", "#FDFCF9", "#CC0606", "#1F78B4"
# (fill, outline): "light" for light backgrounds, "dark" for dark backgrounds
THEMES = {"light": (PAPER, INK), "dark": (INK, PAPER)}
FONT = "'Helvetica Neue', Helvetica, Arial, sans-serif"


def mark(fill: str, line: str) -> str:
    out = []
    for x, y, h, bands in CHROMOSOMES:
        out.append(
            f'<rect x="{x:g}" y="{y:g}" width="{W}" height="{h:g}" rx="{W / 2:g}" fill="{fill}" stroke="{line}" stroke-width="{SW}"/>'
        )
        for by, bh, k in bands:
            out.append(
                f'<rect x="{x + INSET:g}" y="{by:g}" width="{W - 2 * INSET:g}" height="{bh:g}" rx="1.6" fill="{RED if k == "gain" else BLUE}"/>'
            )
        out.append(f'<circle cx="{x + W / 2:g}" cy="{CEN}" r="1.7" fill="{line}"/>')
    return "".join(out)


def icon(fill: str, line: str) -> str:
    vx, vy, vw, vh = VIEW
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vx} {vy} {vw} {vh}" width="{vw * 13}" height="{vh * 13}"><title>cnvscope</title>{mark(fill, line)}</svg>\n'


def wordmark(fill: str, line: str) -> str:
    return (
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="20 22 330 76" width="1320" height="304"><title>cnvscope</title>'
        f'{mark(fill, line)}<text x="112" y="74" font-family="{FONT}" font-weight="500" font-size="44" letter-spacing="-1.3" fill="{line}">cnvscope</text></svg>\n'
    )


def favicon() -> str:
    """Follows the browser's light or dark mode."""
    vx, vy, vw, vh = VIEW
    css = f".c{{fill:{PAPER};stroke:{INK}}}.d{{fill:{INK}}}@media (prefers-color-scheme: dark){{.c{{fill:{INK};stroke:{PAPER}}}.d{{fill:{PAPER}}}}}"
    m = mark("F", "L").replace('fill="F" stroke="L"', 'class="c"').replace('fill="L"', 'class="d"')
    return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vx} {vy} {vw} {vh}"><style>{css}</style>{m}</svg>\n'


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for name, (fill, line) in THEMES.items():
        (OUT / f"cnvscope-icon-{name}.svg").write_text(icon(fill, line))
        (OUT / f"cnvscope-wordmark-{name}.svg").write_text(wordmark(fill, line))
    (OUT.parent / "favicon.svg").write_text(favicon())
    print("wrote", sorted(p.name for p in OUT.iterdir()), "and favicon.svg")
