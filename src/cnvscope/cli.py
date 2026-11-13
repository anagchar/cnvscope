"""Command line.

One command for a run::

    cnvscope --mode dna --tool chisel --input chisel_out/ --outdir results/
    cnvscope --mode dna --chisel --input chisel_out/ --outdir results/      # same

reads the tool's output, writes the dataset and a per-cell table to ``--outdir`` and opens
the viewer. Subcommands for finer control: ``serve``, ``convert``, ``tools``, ``demo``.
"""

from __future__ import annotations

import json
import sys
import threading
import webbrowser
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import click

from . import __version__, readers
from .prepare import prepare

MODES = {"dna": "scDNA", "rna": "scRNA", "atac": "scATAC", "multiome": "multiome"}


@dataclass(frozen=True)
class Tool:
    name: str
    label: str
    mode: str
    #: reader names that can read it; empty = not supported yet
    readers: tuple[str, ...]
    aliases: tuple[str, ...] = ()
    #: other --mode values it serves (Numbat: binned RNA + ATAC multiome runs)
    also: tuple[str, ...] = ()

    @property
    def modes(self) -> tuple[str, ...]:
        return (self.mode, *self.also)


TOOLS = [
    Tool("ascatsc", "ASCAT.sc", "dna", ("ascatsc_json", "ascatsc_r"), ("ascat.sc", "ascat")),
    Tool("hiscanner", "HiScanner", "dna", ("hiscanner",)),
    Tool("chisel", "CHISEL", "dna", ("chisel",)),
    Tool("scphasedcn", "scphasedcn", "dna", ("scphasedcn",)),
    Tool("numbat", "Numbat", "rna", ("numbat",), also=("multiome",)),
    Tool("infercnv", "InferCNV", "rna", ()),
    Tool("copykat", "CopyKAT", "rna", ()),
    Tool("epianeufinder", "epiAneufinder", "atac", ("epianeufinder",), ("epianeu",)),
]
BY_NAME = {k: t for t in TOOLS for k in (t.name, *t.aliases)}


def resolve_tool(name: str) -> Tool:
    t = BY_NAME.get(name.lower())
    if t is None:
        known = ", ".join(t.name for t in TOOLS)
        raise click.UsageError(f"unknown tool {name!r}. Known tools: {known}.")
    return t


def find_input(tool: Tool, path: Path) -> Path:
    """The file or folder the tool's reader takes, from what the user gave as --input."""
    if tool.name == "numbat" and path.is_dir():
        from .readers.numbat import runs_under

        runs = runs_under(path)
        if len(runs) > 1:
            names = ", ".join(str(r.relative_to(path)) for r in runs)
            raise click.UsageError(
                f"{path} holds several Numbat runs: {names}. Point --input at one of them."
            )
        return runs[0] if runs else path
    if tool.name != "ascatsc" or path.is_file():
        return path
    # an ASCAT.sc results folder: prefer the viewer file export_cnvscope() writes, then a web
    # export (both fast), then the R result object (read through R: minutes)
    found = sorted(p for p in path.glob("*.h5ad") if not p.name.startswith("._"))
    found.sort(key=lambda p: not p.name.endswith(".cnvscope.h5ad"))
    if not found:
        found = sorted(path.glob("*.json.gz")) + sorted(path.glob("*.json"))
        found = [p for p in found if readers.REGISTRY["ascatsc_json"].detect(p)]
    if not found:
        objs = [
            p
            for p in path.iterdir()
            if p.suffix.lower() in (".rds", ".rda", ".rdata") and not p.name.startswith("._")
        ]
        objs.sort(key=lambda p: (("result" not in p.name.lower()), -p.stat().st_size))
        found = objs
    if not found:
        raise click.UsageError(f"no ASCAT.sc result (.rds/.Rda or ascat_data.json) in {path}")
    if len(found) > 1:
        click.echo(
            f"cnvscope: using {found[0].name} (also found: {', '.join(p.name for p in found[1:4])})"
        )
    return found[0]


def pick_file(folder: Path | None, file: str) -> tuple[Path, Path]:
    """(input folder, file) for --input-file: a name inside --input, or a path of its own."""
    f = Path(file).expanduser()
    if folder is not None and not f.is_absolute() and (folder / f).exists():
        f = folder / f
    if not f.exists():
        where = f" in {folder}" if folder is not None else ""
        near = ""
        if folder is not None and folder.is_dir():
            names = sorted(p.name for p in folder.iterdir() if not p.name.startswith("."))[:12]
            near = f". Files there: {', '.join(names)}"
        raise click.UsageError(f"--input-file {file} not found{where}{near}")
    return (folder if folder is not None else f.parent), f


def fingerprint(path: Path, depth: int = 2) -> str:
    """Sizes and modification times of the input (a file, or a folder two levels deep)."""
    import hashlib
    import os

    h = hashlib.sha1()

    def walk(d: Path, level: int) -> None:
        with os.scandir(d) as it:
            for e in sorted(it, key=lambda e: e.name):
                if e.name.startswith("."):
                    continue
                if e.is_dir(follow_symlinks=False):
                    if level < depth:
                        walk(Path(e.path), level + 1)
                    continue
                st = e.stat()
                h.update(
                    f"{os.path.relpath(e.path, path)}|{st.st_size}|{st.st_mtime_ns}\n".encode()
                )

    if path.is_file():
        st = path.stat()
        h.update(f"{path.name}|{st.st_size}|{st.st_mtime_ns}".encode())
    else:
        walk(path, 1)
    return h.hexdigest()


def cached_result(outdir: Path, name: str, tool: str, src: Path, fp: str) -> Path | None:
    """The .h5ad of an earlier run on the same, unchanged input, if there is one."""
    run, h5 = outdir / f"{name}_run.json", outdir / f"{name}.h5ad"
    if not (run.exists() and h5.exists()):
        return None
    try:
        meta = json.loads(run.read_text())
    except (OSError, ValueError):
        return None
    same = meta.get("tool") == tool and meta.get("read_from") == str(src.resolve())
    return h5 if same and meta.get("fingerprint") == fp else None


def write_outputs(ds, outdir: Path, name: str, meta: dict) -> dict[str, Path]:
    """<name>.h5ad (open it again with `cnvscope serve`), <name>_cells.tsv, <name>_run.json."""
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "dataset": outdir / f"{name}.h5ad",
        "cells": outdir / f"{name}_cells.tsv",
        "run": outdir / f"{name}_run.json",
    }
    ds.write(paths["dataset"])
    ds.adata.obs.to_csv(paths["cells"], sep="\t", index_label="cell")
    paths["run"].write_text(
        json.dumps(
            meta
            | {
                "cells": int(ds.adata.n_obs),
                "bins": int(ds.adata.n_vars),
                "layers": {k: v.value for k, v in ds.semantics.items()},
                "outputs": {k: str(v) for k, v in paths.items()},
            },
            indent=2,
        )
        + "\n"
    )
    return paths


def add_annotations(ds, path: Path) -> None:
    from .annotate import annotate

    try:
        n = annotate(ds, path)
    except Exception as e:  # noqa: BLE001
        raise click.ClickException(f"could not add annotations from {path}: {e}") from e
    click.echo(f"cnvscope: annotations for {n:,} of {ds.adata.n_obs:,} cells from {path.name}")


def free_port(host: str, port: int, tries: int = 50) -> int:
    """`port`, or the next free one after it: another cnvscope may already be running."""
    import socket

    for p in range(port, port + tries):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind((host, p))
            except OSError:
                continue
            return p
    raise click.ClickException(f"no free port between {port} and {port + tries - 1}; use --port")


def run_viewer(ds, name: str, host: str, port: int, no_open: bool) -> None:
    import uvicorn

    from .server import create_app

    chosen = free_port(host, port)
    if chosen != port:
        click.echo(
            f"cnvscope: port {port} is in use (another viewer still running?), using {chosen}"
        )
    port = chosen
    url = f"http://{host}:{port}"
    click.echo(f"cnvscope: viewer at {url}  (Ctrl+C to stop)")
    if not no_open:
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    uvicorn.run(create_app(ds, name), host=host, port=port, log_level="warning")


def tool_flag(t: Tool):
    names = [f"--{t.name}"] + [f"--{a}" for a in t.aliases if "." not in a]
    return click.option(*names, f"flag_{t.name}", is_flag=True, help=f"Same as --tool {t.name}.")


def with_tool_flags(f):
    for t in reversed(TOOLS):
        f = tool_flag(t)(f)
    return f


@click.group(invoke_without_command=True, context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--mode", type=click.Choice(list(MODES)), help="Data type: dna, rna, atac or multiome."
)
@click.option(
    "--tool",
    "tool_name",
    metavar="NAME",
    help="Tool that made --input: " + ", ".join(t.name for t in TOOLS) + ".",
)
@with_tool_flags
@click.option(
    "--input",
    "-i",
    "input_path",
    type=click.Path(exists=True, path_type=Path),
    help="The tool's output folder (or file).",
)
@click.option(
    "--input-file",
    "--input-file-name",
    "--input_file",
    "input_file",
    metavar="FILE",
    help="A specific file to read: a name inside --input or a full path, e.g. "
    "--input-file results.rds for an ASCAT.sc folder holding several result objects.",
)
@click.option(
    "--outdir",
    "-o",
    type=click.Path(path_type=Path),
    help="Folder for the dataset (.h5ad), the per-cell table and the run summary.",
)
@click.option(
    "--annotations",
    "-a",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Cell annotations (.csv/.tsv): a barcode column (cell/barcode, or the first) and any "
    "others, e.g. cell type; shown as tracks beside the heatmap.",
)
@click.option("--name", help="Name of the dataset and its files (default: the input's name).")
@click.option("--no-view", is_flag=True, help="Only write the outputs; do not start the viewer.")
@click.option(
    "--force",
    is_flag=True,
    help="Read the input again even if --outdir already has a result for it.",
)
@click.option("--port", default=8765, show_default=True, help="Viewer port.")
@click.option("--host", default="127.0.0.1", show_default=True, help="Viewer address.")
@click.option("--no-open", is_flag=True, help="Do not open a browser.")
@click.version_option(__version__, "-V", "--version", prog_name="cnvscope")
@click.pass_context
def main(
    ctx: click.Context,
    mode,
    tool_name,
    input_path,
    input_file,
    outdir,
    annotations,
    name,
    no_view,
    force,
    port,
    host,
    no_open,
    **flags,
) -> None:
    """cnvscope: interactive viewer for single-cell copy-number calls.

    \b
    One command:
      cnvscope --mode dna --tool chisel --input chisel_out/ --outdir results/
      cnvscope --mode dna --hiscanner --input hiscanner_out/ --outdir results/
      cnvscope --mode dna --ascatsc --input ascatsc_results/ --outdir results/
      cnvscope --mode atac --epianeufinder --input epi_out/ --annotations cells.csv \\
               --outdir results/
      cnvscope --mode multiome --numbat --input numbat_run/ --outdir results/

    \b
    Open the outputs again later:
      cnvscope serve results/<name>.h5ad
    """
    if ctx.invoked_subcommand is not None:
        return
    chosen = [n.removeprefix("flag_") for n, on in flags.items() if on]
    if tool_name:
        chosen.append(tool_name)
    named_by = input_path  # the folder names the dataset; a lone --input-file names it itself
    if input_file:
        input_path, file = pick_file(input_path, input_file)
    if not (mode or chosen or input_path or outdir):
        click.echo(ctx.get_help())
        return
    missing = [
        o
        for o, v in (
            ("--mode", mode),
            ("--tool", chosen),
            ("--input (or --input-file with a full path)", input_path),
            ("--outdir", outdir),
        )
        if not v
    ]
    if missing:
        raise click.UsageError(
            f"missing {', '.join(missing)}. Example: "
            "cnvscope --mode dna --tool chisel --input chisel_out/ --outdir results/"
        )
    tools = {resolve_tool(c).name for c in chosen}
    if len(tools) > 1:
        raise click.UsageError(f"choose one tool, got {', '.join(sorted(tools))}")
    tool = resolve_tool(chosen[0])
    if mode not in tool.modes:
        fits = ", ".join(t.name for t in TOOLS if mode in t.modes) or "none yet"
        raise click.UsageError(
            f"{tool.label} is a {' / '.join(MODES[m] for m in tool.modes)} tool, not "
            f"{MODES[mode]}. Tools for --mode {mode}: {fits}."
        )
    if not tool.readers:
        ready = ", ".join(t.name for t in TOOLS if t.readers)
        raise click.UsageError(
            f"{tool.label} is not supported yet (planned). Supported now: {ready}."
        )

    src = file if input_file else find_input(tool, input_path)
    opts = {}
    if src.suffix == ".h5ad":  # a viewer file written earlier (export_cnvscope or convert)
        from .model import CNVDataset

        reader = readers.Reader("h5ad", MODES[mode], lambda p: True, CNVDataset.read)
    else:
        reader = next(
            (readers.REGISTRY[r] for r in tool.readers if readers.REGISTRY[r].detect(src)), None
        )
        if reader is not None and reader.name == "ascatsc_r":
            opts["use_viewer_file"] = not force  # --force: read the R object itself
        if reader is not None and reader.name == "numbat":
            opts["modality"] = MODES[mode]
    if reader is None:
        raise click.UsageError(
            f"{src} does not look like {tool.label} output. "
            "Point --input at the tool's output folder."
        )
    name = name or (named_by or src).resolve().name.split(".")[0]
    fp = fingerprint(src)
    cached = None if force else cached_result(outdir, name, tool.name, src, fp)
    if cached is not None:
        # same input, unchanged since the last run: reopen its result (seconds, not minutes)
        click.echo(f"cnvscope: {src.name} is unchanged since the last run; reopening {cached}")
        click.echo("          (use --force to read it again)")
        from .model import CNVDataset

        ds = prepare(CNVDataset.read(cached))
        if annotations:
            add_annotations(ds, annotations)
            ds.write(cached)  # keep the outputs in step with the annotations
            ds.adata.obs.to_csv(outdir / f"{name}_cells.tsv", sep="\t", index_label="cell")
        click.echo(f"cnvscope: {ds.adata.n_obs:,} cells × {ds.adata.n_vars:,} bins")
        if not no_view:
            run_viewer(ds, name, host, port, no_open)
        return
    click.echo(f"cnvscope: reading {tool.label} output from {src}")
    try:
        ds = prepare(reader.read(src, **opts))
    except Exception as e:  # noqa: BLE001 - report reader failures plainly
        raise click.ClickException(f"could not read {src}: {e}") from e
    if annotations:
        add_annotations(ds, annotations)
    meta = {
        "cnvscope": __version__,
        "mode": mode,
        "tool": tool.name,
        "input": str(input_path.resolve()),
        "read_from": str(src.resolve()),
        "fingerprint": fp,
        "annotations": str(annotations.resolve()) if annotations else None,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command": " ".join(sys.argv),
    }
    paths = write_outputs(ds, outdir, name, meta)
    click.echo(f"cnvscope: {ds.adata.n_obs:,} cells × {ds.adata.n_vars:,} bins")
    for k, p in paths.items():
        click.echo(f"  {k:8} {p}")
    if not no_view:
        run_viewer(ds, name, host, port, no_open)


@main.command()
def rscript() -> None:
    """Print the path of export_cnvscope.R, for the end of an ASCAT.sc pipeline.

    \b
    In R:
      source(system("cnvscope rscript", intern = TRUE))
      export_cnvscope(res, "results/result_object.Rda")
    writes results/result_object.cnvscope.h5ad, which cnvscope then opens instead of the .Rda.
    """
    click.echo(Path(__file__).parent / "R" / "export_cnvscope.R")


@main.command()
def tools() -> None:
    """List the tools and which ones can be read now."""
    for t in TOOLS:
        state = "ready" if t.readers else "planned"
        click.echo(f"  --mode {'/'.join(t.modes):14} --tool {t.name:14} {t.label:14} {state}")


@main.command()
@click.option("-o", "--output", type=click.Path(path_type=Path), default="demo.h5ad")
@click.option("--cells", default=300, show_default=True)
@click.option("--seed", default=0, show_default=True)
def demo(output: Path, cells: int, seed: int) -> None:
    """Write a synthetic scDNA dataset."""
    from .readers.demo import make_demo

    make_demo(n_cells=cells, seed=seed).write(output)
    click.echo(f"wrote {output}")


@main.command()
@click.argument("path", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", type=click.Path(path_type=Path), required=True)
@click.option("--tool", help="Reader to use (detected if omitted).")
@click.option("--no-cluster", is_flag=True, help="Keep the tool's cell order.")
def convert(path: Path, output: Path, tool: str | None, no_cluster: bool) -> None:
    """Convert a tool's output to a cnvscope .h5ad."""
    ds = prepare(readers.load(path, tool), cluster=not no_cluster)
    ds.write(output)
    click.echo(f"wrote {output}: {ds.adata.n_obs} cells × {ds.adata.n_vars} bins ({ds.tool})")


@main.command()
@click.argument("path", type=click.Path(path_type=Path), required=False)
@click.option("--tool", help="Reader to use (detected if omitted).")
@click.option("--host", default="127.0.0.1", show_default=True)
@click.option("--port", default=8765, show_default=True)
@click.option("--no-open", is_flag=True, help="Do not open a browser.")
def serve(path: Path | None, tool: str | None, host: str, port: int, no_open: bool) -> None:
    """Open PATH (.h5ad, a tool's output, or `demo`) in the viewer.

    Without PATH the viewer starts on its upload screen.
    """
    ds, name = None, ""
    if path is not None:
        if str(path) != "demo" and not path.exists():
            raise click.BadParameter(f"{path} does not exist", param_hint="PATH")
        ds = prepare(readers.load(path, tool))
        name = "Demo dataset" if str(path) == "demo" else path.name.split(".")[0]
        click.echo(f"cnvscope: {ds.adata.n_obs} cells × {ds.adata.n_vars} bins ({ds.tool})")
    run_viewer(ds, name, host, port, no_open)


if __name__ == "__main__":
    main()
