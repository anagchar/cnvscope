"""scphasedcn caller folders (normal and tumour mode) and the extraction's 100 kb bins."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from click.testing import CliRunner

from cnvscope.cli import main
from cnvscope.readers import detect, load

CELLS = ["c1", "c2", "c3"]
REGIONS = pd.DataFrame(
    {
        "region": ["chr1", "chr2", "chrX", "chr1p"],
        "class": ["chrom", "chrom", "chrom", "arm"],
        "chrom": ["chr1", "chr2", "chrX", "chr1"],
        "start": [0, 0, 0, 0],
        "end": [6e6, 4e6, 3e6, 2.4e6],
    }
)


def caller(root: Path, sex: str = "female") -> Path:
    d = root / "caller"
    d.mkdir(parents=True)
    pd.DataFrame({"cell": CELLS, "qc": True, "cell_type": ["a", "b", "a"], "fold": np.nan}).to_csv(
        d / "cells.tsv", sep="\t", index=False
    )
    REGIONS.to_csv(d / "regions.tsv", sep="\t", index=False)
    (d / "thresholds.tsv").write_text("score\tclass\tthr\n")
    (d / "params.yaml").write_text("mode: normal\n")
    (d / "NOTES.md").write_text(f"# run\n\n- dir: extract\n- sex: {sex}\n")
    return d


def normal_run(root: Path, sex: str = "female") -> Path:
    d = caller(root, sex)
    calls = pd.DataFrame(
        # c1: chr1 gain, but a stronger nested chr1p loss; c2: chrX loss; a pool row and an
        # unreportable row that must not be painted
        {
            "cell": ["c1", "c1", "c2", "pool1", "c3"],
            "region": ["chr1", "chr1p", "chrX", "chr2", "chr2"],
            "class": ["chrom", "arm", "chrom", "chrom", "chrom"],
            "chrom": ["chr1", "chr1", "chrX", "chr2", "chr2"],
            "start": [0, 0, 0, 0, 0],
            "end": [6e6, 2.4e6, 3e6, 4e6, 4e6],
            "state": ["gain", "loss", "loss", "gain", "cnloh"],
            "called": True,
            "level": ["cell", "cell", "cell", "pool", "cell"],
            "source": "cell",
            "reportable": [True, True, True, True, False],
        }
    )
    calls.to_csv(d / "calls.tsv", sep="\t", index=False)
    return root


def tumour_run(root: Path) -> Path:
    d = caller(root, "male")
    cells = pd.read_csv(d / "cells.tsv", sep="\t")
    cells["clone"] = ["C1", "N0", "C1"]
    cells["assigned"] = [True, True, False]
    cells.to_csv(d / "cells.tsv", sep="\t", index=False)
    (d / "calls.tsv").write_text("cell\tregion\tclass\tchrom\tstart\tend\tstate\tcalled\tlevel\n")
    seg = pd.DataFrame(
        {
            "clone": ["C1", "C1", "C1", "C1", "N0", "N0", "N0"],
            "chrom": ["chr1", "chr1", "chr2", "chrX", "chr1", "chr2", "chrX"],
            "start": [0, 3e6, 0, 0, 0, 0, 0],
            "end": [3e6, 6e6, 4e6, 3e6, 6e6, 4e6, 3e6],
            "state": ["3:1", "1:1", "2:0", "1:0", "1:1", "1:1", "1:0"],
        }
    )
    seg.to_csv(d / "clone_segments.tsv", sep="\t", index=False)
    return root


def extract(root: Path) -> None:
    """bins_100kb.h5ad: 100 kb bins, sparse counts and haplotype votes."""
    d = root / "extract"
    d.mkdir(parents=True)
    var = pd.concat(
        pd.DataFrame({"chrom": c, "start": np.arange(0, n, 1e5), "end": np.arange(1e5, n + 1, 1e5)})
        for c, n in (("chr1", 6e6), ("chr2", 4e6), ("chrX", 3e6))
    ).reset_index(drop=True)
    var.index = var.index.astype(str)
    nb = len(var)
    X = np.ones((4, nb))  # an extra cell the caller left out
    hapA, hapB = np.full((4, nb), 3.0), np.ones((4, nb))
    a = ad.AnnData(
        X=sp.csr_matrix(X),
        obs=pd.DataFrame(index=["c3", "c1", "c2", "c0"]),
        var=var,
        layers={"hapA": sp.csr_matrix(hapA), "hapB": sp.csr_matrix(hapB)},
    )
    a.write_h5ad(d / "bins_100kb.h5ad")


def win(ds, chrom: str, mb: float) -> int:
    v = ds.adata.var
    return int(
        np.flatnonzero((v["chr"] == chrom) & (v["start"] <= mb * 1e6) & (v["end"] > mb * 1e6))[0]
    )


def test_normal_mode_paints_reportable_calls(tmp_path):
    root = normal_run(tmp_path / "s")
    assert detect(root).name == "scphasedcn" and detect(root / "caller").name == "scphasedcn"
    ds = load(root)
    a = ds.adata
    assert a.uns["scphasedcn_mode"] == "normal" and a.n_vars == 6 + 4 + 3
    assert a.obs_names.tolist() == CELLS and "fold" not in a.obs  # empty columns dropped
    T, M, m = ds.layer("total_cn"), ds.layer("major_cn"), ds.layer("minor_cn")
    # c1: the nested 1p loss (1:0) over the chromosome gain (2:1)
    assert (T[0, win(ds, "chr1", 1)], M[0, win(ds, "chr1", 1)], m[0, win(ds, "chr1", 1)]) == (
        1,
        1,
        0,
    )
    assert T[0, win(ds, "chr1", 4)] == 3
    assert T[1, win(ds, "chrX", 1)] == 1  # c2: female chrX loss
    assert (T[2] == 2).all()  # the unreportable row and the pool row are not painted


def test_male_x_is_haploid(tmp_path):
    ds = load(normal_run(tmp_path / "s", sex="male"))
    x = win(ds, "chrX", 1)
    T = ds.layer("total_cn")
    assert T[0, x] == 1 and T[2, x] == 1  # neutral = one copy
    assert T[1, x] == 0  # a loss of the only copy


def test_tumour_mode_paints_clone_segments(tmp_path):
    ds = load(tumour_run(tmp_path / "t"))
    T, M, m = ds.layer("total_cn"), ds.layer("major_cn"), ds.layer("minor_cn")
    assert ds.adata.uns["scphasedcn_mode"] == "tumour"
    assert (T[0, win(ds, "chr1", 1)], M[0, win(ds, "chr1", 1)], m[0, win(ds, "chr1", 1)]) == (
        4,
        3,
        1,
    )
    assert (M[0, win(ds, "chr2", 1)], m[0, win(ds, "chr2", 1)]) == (2, 0)
    assert T[1, win(ds, "chr1", 1)] == 2  # N0
    assert np.isnan(T[2]).all()  # not assigned
    assert ds.adata.obs["clone"].tolist() == ["C1", "N0", "C1"]


def test_extraction_adds_read_counts_and_haplotype_fraction(tmp_path):
    root = normal_run(tmp_path / "s")
    extract(root)
    ds = load(root)
    R = ds.layer("read_counts")
    assert (R == 10).all()  # ten 100 kb bins per 1 Mb window, cells matched by name
    np.testing.assert_allclose(ds.layer("hap_a_frac"), 0.75)
    assert ds.semantics["hap_a_frac"].value == "baf"


def test_cli_run(tmp_path):
    root = tumour_run(tmp_path / "t")
    r = CliRunner().invoke(
        main,
        ["--mode", "dna", "--scphasedcn", "-i", str(root), "-o", str(tmp_path / "o"), "--no-view"],
    )
    assert r.exit_code == 0, r.output
    a = ad.read_h5ad(tmp_path / "o" / "t.h5ad")
    assert a.uns["tool"] == "scphasedcn" and set(a.layers.keys()) >= {"total_cn", "major_cn"}


def _r_with_data_table() -> bool:
    import shutil
    import subprocess

    exe = shutil.which("Rscript")
    if not exe:
        return False
    r = subprocess.run([exe, "-e", "library(data.table)"], capture_output=True)
    return r.returncode == 0


@pytest.mark.skipif(not _r_with_data_table(), reason="needs R with data.table")
def test_result_object_is_read_through_r(tmp_path):
    """<sample>_call_scn.rds (call_scn() saved with saveRDS) beside <sample>/extract."""
    import subprocess

    src = normal_run(tmp_path / "tables")
    d = src / "caller"
    root = tmp_path / "scphasedcn"
    (root / "S1").mkdir(parents=True)
    extract(root / "S1")
    rds = root / "S1_call_scn.rds"
    subprocess.run(
        [
            "Rscript",
            "-e",
            "library(data.table); a <- commandArgs(TRUE); d <- a[1]; "
            "x <- list(mode = 'normal', cells = fread(file.path(d, 'cells.tsv')), "
            "calls = fread(file.path(d, 'calls.tsv')), "
            "regions = fread(file.path(d, 'regions.tsv')), "
            "thresholds = data.table(), run = list(inputs = c('dir: /not/here', 'sex: female'))); "
            "class(x) <- 'scphasedcn'; saveRDS(x, a[2])",
            str(d),
            str(rds),
        ],
        check=True,
        capture_output=True,
    )
    for p in (rds, root / "S1"):  # the object itself, or the sample folder beside it
        ds = load(p, "scphasedcn")
        T = ds.layer("total_cn")
        assert T[0, win(ds, "chr1", 1)] == 1 and T[0, win(ds, "chr1", 4)] == 3
        assert (ds.layer("read_counts") == 10).all()  # S1/extract found beside the object
        assert ds.adata.uns["source"] == "S1"


def test_extract_alone_is_explained(tmp_path):
    extract(tmp_path / "S1")
    with pytest.raises(ValueError, match="extract/ holds only the read counts"):
        load(tmp_path / "S1", "scphasedcn")
