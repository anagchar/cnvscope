"""Tool readers: each turns one tool's output into a :class:`CNVDataset`.

A reader module calls :func:`register` with a name, a ``detect(path) -> bool`` that recognises
the tool's output (a file or a directory) cheaply, and a ``read(path, **opts) -> CNVDataset``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from ..model import CNVDataset


@dataclass(frozen=True)
class Reader:
    name: str
    modality: str
    detect: Callable[[Path], bool]
    read: Callable[..., CNVDataset]
    description: str = ""


REGISTRY: dict[str, Reader] = {}


def register(reader: Reader) -> Reader:
    REGISTRY[reader.name] = reader
    return reader


def detect(path: str | Path) -> Reader | None:
    p = Path(path)
    for r in REGISTRY.values():
        if r.name == "demo":
            continue
        try:
            if r.detect(p):
                return r
        except OSError:
            continue
    return None


def load(path: str | Path, tool: str | None = None, **opts) -> CNVDataset:
    """Read ``path``: an .h5ad written by cnvscope, or any registered tool's output."""
    p = Path(path)
    if p.suffix == ".h5ad":
        return CNVDataset.read(p)
    if str(p) == "demo" and not p.exists():
        return REGISTRY["demo"].read(p, **opts)
    if tool:
        if tool not in REGISTRY:
            raise ValueError(f"unknown tool {tool!r}; known: {', '.join(sorted(REGISTRY))}")
        return REGISTRY[tool].read(p, **opts)
    r = detect(p)
    if r is None:
        raise ValueError(
            f"could not recognise {p} as the output of any of: {', '.join(sorted(REGISTRY))}. "
            "Pass --tool to choose one."
        )
    return r.read(p, **opts)


# Import reader modules so they register themselves.
from . import (  # noqa: E402,F401
    ascatsc_json,
    ascatsc_r,
    chisel,
    demo,
    epianeufinder,
    hiscanner,
    numbat,
    scphasedcn,
)
