"""epiAneufinder's results_table.tsv, cell annotations, and the --mode atac run."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from cnvscope.annotate import annotate
from cnvscope.cli import main
from cnvscope.model import Semantics
from cnvscope.readers import detect, load

CELLS = ["AAACGAAAGC-1", "TTGCAGGTCA-1", "GGGTCATTAC-1"]
# chr1 windows with one filtered out (3-4 Mb), and chr2
WINDOWS = [("chr1", 1, 1e6), ("chr1", 1e6 + 1, 2e6), ("chr1", 2e6 + 1, 3e6), ("chr1", 4e6 + 1, 5e6)]
WINDOWS += [("chr2", 1, 1e6), ("chr2", 1e6 + 1, 2e6)]
CALLS = np.array(  # windows x cells, 0 loss / 1 normal / 2 gain
    [[1, 2, 0], [1, 2, 0], [1, 2, 1], [1, 1, 1], [1, 1, 0], [1, 1, 0]]
)


def r_num(v: float) -> str:
    # R's write.table prints 1e+06 for one million
    return f"{v:.0e}" if v >= 1e5 and v % 1e5 == 0 else str(int(v))


def write_epi(root: Path) -> Path:
    """An epiAneufinder output folder, written like write.table(..., quote = FALSE)."""
    d = root / "epiAneufinder_results"
    d.mkdir(parents=True)
    head = "seq start end " + " ".join(f"cell-{c}" for c in CELLS)
    rows = [
        f"{i + 1} {c} {r_num(s)} {r_num(e)} " + " ".join(map(str, CALLS[i]))
        for i, (c, s, e) in enumerate(WINDOWS)
    ]
    (d / "results_table.tsv").write_text("\n".join([head, *rows]) + "\n")
    return root


def test_reads_results_table(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    assert detect(root).name == "epianeufinder"
    ds = load(root)
    a = ds.adata
    assert a.uns["tool"] == "epiAneufinder" and a.uns["modality"] == "scATAC"
    assert a.obs_names.tolist() == CELLS
    assert ds.semantics == {"state": Semantics.TERNARY_STATE}
    assert a.var["end"].tolist()[:2] == [1_000_000, 2_000_000]
    np.testing.assert_array_equal(ds.layer("state"), CALLS.T - 1)
    np.testing.assert_allclose(a.obs["cnv_burden"], [0, 3 / 6, 4 / 6])
    np.testing.assert_allclose(a.obs["gain_fraction"], [0, 3 / 6, 0])
    np.testing.assert_allclose(a.obs["loss_fraction"], [0, 0, 4 / 6])


def test_reads_the_table_file_and_results_folder(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    for p in (root / "epiAneufinder_results", root / "epiAneufinder_results" / "results_table.tsv"):
        assert load(p, "epianeufinder").adata.n_obs == 3


@pytest.mark.parametrize(
    "barcodes",
    [CELLS, [f"cell-{c}" for c in CELLS], [f"sample1#{c[:-2]}" for c in CELLS]],
)
def test_annotations_match_barcodes(tmp_path, barcodes):
    ds = load(write_epi(tmp_path / "epi"))
    f = tmp_path / "annot.csv"
    pd.DataFrame({"cell": barcodes, "annot": ["normal", "tumour", "tumour"]}).to_csv(f, index=False)
    assert annotate(ds, f) == 3
    assert ds.adata.obs["cell_type"].tolist() == ["normal", "tumour", "tumour"]


def test_annotations_that_match_nothing_are_an_error(tmp_path):
    ds = load(write_epi(tmp_path / "epi"))
    f = tmp_path / "annot.tsv"
    f.write_text("barcode\tcelltype\nXXXX\tB cell\n")
    with pytest.raises(ValueError, match="no barcode"):
        annotate(ds, f)


def test_atac_run_with_annotations(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    f = tmp_path / "cells.tsv"
    f.write_text(
        "barcode\tcell_type\tsample\n" + "".join(f"{c}\tT{i}\tS1\n" for i, c in enumerate(CELLS))
    )
    out = tmp_path / "results"
    r = CliRunner().invoke(
        main,
        [
            "--mode",
            "atac",
            "--epianeufinder",
            "--input",
            str(root),
            "--annotations",
            str(f),
            "--outdir",
            str(out),
            "--no-view",
        ],
    )
    assert r.exit_code == 0, r.output
    assert "annotations for 3 of 3 cells" in r.output
    a = ad.read_h5ad(out / "epi_out.h5ad")
    assert a.obs["cell_type"].tolist() == ["T0", "T1", "T2"] and a.uns["tool"] == "epiAneufinder"
    cells = pd.read_csv(out / "epi_out_cells.tsv", sep="\t")
    assert {"cnv_burden", "cell_type", "sample"} <= set(cells.columns)


def test_atac_tool_needs_atac_mode(tmp_path):
    root = write_epi(tmp_path / "epi_out")
    r = CliRunner().invoke(
        main,
        ["--mode", "dna", "--tool", "epianeufinder", "-i", str(root), "-o", str(tmp_path / "o")],
    )
    assert r.exit_code != 0 and "scATAC tool" in r.output
