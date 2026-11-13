"""Numbat out_dir (binned multiome run): final iteration, states at p >= 0.9, CLI modes."""

import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from cnvscope.cli import main
from cnvscope.readers import detect, load

CELLS = ["AAAC", "GGTT", "TTAG", "CCGA"]
# 1 Mb bins on chr1 (0-5 Mb) and chr2 (0-3 Mb), as Numbat-multiome names them
BINS = [("chr1", k * 1_000_000, (k + 1) * 1_000_000) for k in range(5)]
BINS += [("chr2", k * 1_000_000, (k + 1) * 1_000_000) for k in range(3)]


def write_run(d: Path, it: int = 2) -> Path:
    d.mkdir(parents=True)
    for i in range(1, it + 1):
        pd.DataFrame(
            {
                "cell": CELLS,
                "clone_opt": [1, 2, 2, 1],
                "GT_opt": ["", "1a", "1a,2a", ""],
                "p_opt": [1.0, 0.99, 0.98, 1.0],
                "p_cnv": [1e-9, 1.0, 1.0, 0.2],
                "compartment_opt": ["normal", "tumor", "tumor", "normal"],
            }
        ).to_csv(d / f"clone_post_{i}.tsv", sep="\t", index=False)
        pd.DataFrame(
            {
                "sample": 2,
                "CHROM": [1, 1, 2],
                "seg": ["1b", "1d", "2b"],
                "cnv_state": ["del", "amp", "loh"],
                "seg_start": [500_000, 1_200_000, 100_000],  # 1a overlaps bins 0-1; 1b bins 1-2
                "seg_end": [1_600_000, 2_900_000, 2_500_000],
                "seg_cons": ["1a", "1b", "2a"],
            }
        ).to_csv(d / f"segs_consensus_{i}.tsv", sep="\t", index=False)
        pd.DataFrame(
            {
                "cell": CELLS[:3],  # CCGA has no genotype
                "1a": [0.01, 0.95, 0.60],
                "1b": [0.02, 0.97, 0.99],
                "2a": [0.50, 0.10, 0.92],
            }
        ).to_csv(d / f"geno_{i}.tsv", sep="\t", index=False)
    names = [f"{c}:{s}-{e}" for c, s, e in BINS]
    rng = np.random.default_rng(0)
    g = pd.DataFrame(rng.normal(0, 0.1, (4, len(BINS))), columns=names)
    g.insert(0, "cell", CELLS)
    with gzip.open(d / "gexp_roll_wide.tsv.gz", "wt") as f:
        g[["cell", *names[::-1]]].to_csv(f, sep="\t", index=False)  # any column order
    return d


def test_reads_final_iteration(tmp_path):
    d = write_run(tmp_path / "run" / "numbat_out")
    for p in (d, d.parent):
        assert detect(p).name == "numbat"
    ds = load(d.parent)
    a = ds.adata
    assert a.uns["tool"] == "Numbat" and a.uns["numbat_iteration"] == 2
    assert a.uns["modality"] == "multiome" and a.uns["numbat_binned"]
    assert a.n_vars == 8 and a.var["chr"].tolist()[:2] == ["chr1", "chr1"]
    assert a.obs["clone"].tolist() == ["1", "2", "2", "1"]
    assert a.obs["compartment"].tolist() == ["normal", "tumor", "tumor", "normal"]
    assert set(ds.layer_names) == {"total_cn", "major_cn", "minor_cn", "p_cnv", "logr"}


def test_states_painted_where_numbat_is_confident(tmp_path):
    ds = load(write_run(tmp_path / "numbat_out"))
    P, T = ds.layer("p_cnv"), ds.layer("total_cn")
    M, m = ds.layer("major_cn"), ds.layer("minor_cn")
    # GGTT: 1a del (0.95) on bin 0; bin 1 is under 1a and 1b: the more likely one (1b amp) wins
    assert P[1, 0] == pytest.approx(0.95) and (M[1, 0], m[1, 0]) == (1, 0)
    assert P[1, 1] == pytest.approx(0.97) and T[1, 1] == 3
    # TTAG: 1a at 0.60 is below 0.9 -> neutral; 2a loh at 0.92 -> 2:0
    assert P[2, 0] == pytest.approx(0.60) and T[2, 0] == 2 and m[2, 0] == 1
    assert (M[2, 5], m[2, 5]) == (2, 0)
    assert P[2, 4] == 0  # chr1 4-5 Mb: no CNV segment
    assert (T[0] == 2).all()  # AAAC: nothing confident
    assert np.isnan(T[3]).all() and np.isnan(P[3]).all()  # no genotype
    # TTAG: the chr1 amp (2 bins) and the chr2 loh (3 bins)
    np.testing.assert_allclose(ds.adata.obs["cnv_burden"][:3], [0, 3 / 8, 5 / 8])


def test_expression_follows_the_bins(tmp_path):
    d = write_run(tmp_path / "numbat_out")
    g = pd.read_csv(d / "gexp_roll_wide.tsv.gz", sep="\t").set_index("cell")
    ds = load(d)
    np.testing.assert_allclose(ds.layer("logr")[1, 0], g.loc["GGTT", "chr1:0-1000000"], rtol=1e-6)


def test_no_dna_qc_metrics(tmp_path):
    from cnvscope.prepare import prepare

    ds = prepare(load(write_run(tmp_path / "numbat_out")))
    assert "mapd" not in ds.adata.obs and "median_residual" not in ds.adata.obs
    assert "total" in ds.dendrograms  # clustered on the posteriors


def run(*args):
    return CliRunner().invoke(main, [*args, "--no-view"])


def test_cli_modes(tmp_path):
    d = write_run(tmp_path / "s" / "numbat_out")
    r = run("--mode", "multiome", "--numbat", "-i", str(d.parent), "-o", str(tmp_path / "o"))
    assert r.exit_code == 0, r.output
    assert ad.read_h5ad(tmp_path / "o" / "s.h5ad").uns["modality"] == "multiome"
    r = run("--mode", "rna", "--tool", "numbat", "-i", str(d), "-o", str(tmp_path / "o2"))
    assert r.exit_code == 0, r.output
    assert ad.read_h5ad(tmp_path / "o2" / "numbat_out.h5ad").uns["modality"] == "scRNA"
    r = run("--mode", "dna", "--numbat", "-i", str(d), "-o", str(tmp_path / "o3"))
    assert r.exit_code != 0 and "scRNA / multiome tool" in r.output


def test_cli_several_runs(tmp_path):
    write_run(tmp_path / "all" / "1MB" / "numbat_out")
    write_run(tmp_path / "all" / "220KB" / "numbat_out")
    r = run(
        "--mode", "multiome", "--numbat", "-i", str(tmp_path / "all"), "-o", str(tmp_path / "o")
    )
    assert r.exit_code != 0 and "several Numbat runs" in r.output and "220KB" in r.output
