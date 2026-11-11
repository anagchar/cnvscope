"""ASCAT.sc result objects from R (``.rds``, or ``.Rda``/``.RData`` holding one).

Converted with the ASCAT.sc Explorer's own R exporter (``R/ascatsc_to_web.R``, shipped in this
package) through ``Rscript``, then read as an ASCAT.sc web export. Needs R with the
``data.table`` and ``jsonlite`` packages. Large results take minutes, so a viewer file beside
the object is used instead when there is one (see :func:`viewer_file`): write it once at the end
of the ASCAT.sc pipeline with ``export_cnvscope()`` (``R/export_cnvscope.R``).
"""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from ..model import CNVDataset
from . import Reader, register
from .ascatsc_json import _open, from_dict

SCRIPT = Path(__file__).parents[1] / "R" / "export_ascatsc.R"
SUFFIXES = (".rds", ".rda", ".rdata")
#: viewer files looked for beside results/x.Rda, best first
VIEWER_SUFFIXES = (".cnvscope.h5ad", ".h5ad", ".cnvscope.json.gz", ".cnvscope.json")


def viewer_file(path: Path) -> Path | None:
    """The small file cnvscope opens instead of the R object ``path``, if one is beside it.

    ``results/x.Rda`` -> ``results/x.cnvscope.h5ad`` (written by ``export_cnvscope()``),
    ``results/x.h5ad``, or ``results/x.cnvscope.json.gz``. One older than the object is
    ignored: the object was saved again after it.
    """
    stem = path.name[: -len(path.suffix)] if path.suffix.lower() in SUFFIXES else path.name
    try:
        made = path.stat().st_mtime
    except OSError:
        return None
    for suf in VIEWER_SUFFIXES:
        f = path.with_name(stem + suf)
        try:
            if f.is_file() and f.stat().st_mtime >= made - 1:
                return f
            if f.is_file():
                print(f"cnvscope: ignoring {f.name}, it is older than {path.name}", file=sys.stderr)
        except OSError:
            continue
    return None


def read_viewer_file(f: Path) -> CNVDataset:
    if f.suffix == ".h5ad":
        return CNVDataset.read(f)
    return from_dict(_open(f), name=f.name.split(".")[0])


def _detect(path: Path) -> bool:
    return path.is_file() and path.suffix.lower() in SUFFIXES


def read_ascatsc_r(
    path: Path, rscript: str | None = None, use_viewer_file: bool = True
) -> CNVDataset:
    """Read an ASCAT.sc .rds/.Rda: from its viewer file if it has one, else through R."""
    if use_viewer_file and (f := viewer_file(path)) is not None:
        print(f"cnvscope: using {f.name} beside {path.name} (fast; no R needed)", file=sys.stderr)
        return read_viewer_file(f)
    exe = rscript or shutil.which("Rscript")
    if not exe:
        raise RuntimeError(
            "reading an ASCAT.sc .rds/.Rda needs R (Rscript on the PATH). Or export it in R with "
            "rds_to_web() and open the .json."
        )
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "ascat_data.json"
        r = subprocess.run([exe, str(SCRIPT), str(path), str(out)], capture_output=True, text=True)
        if r.returncode != 0 or not out.exists():
            tail = "\n".join((r.stderr or r.stdout).strip().splitlines()[-15:])
            raise RuntimeError(f"R export of {path.name} failed:\n{tail}")
        d = _open(out)
    return from_dict(d, name=path.name.rsplit(".", 1)[0])


def _read(path: Path, **opts) -> CNVDataset:
    return read_ascatsc_r(path, **opts)


register(
    Reader("ascatsc_r", "scDNA", _detect, _read, "ASCAT.sc result object (.rds / .Rda), via R")
)
