"""An ASCAT.sc .rds/.Rda is opened from the viewer file beside it, without R."""

import os
import shutil
from pathlib import Path

import pytest
from click.testing import CliRunner

from cnvscope.cli import main
from cnvscope.readers import ascatsc_r, load

FIX = Path(__file__).parents[1] / "fixtures" / "ascatsc_json" / "demo_export.json.gz"


@pytest.fixture
def no_r(monkeypatch):
    # any attempt to go through R fails, so a pass proves R was not needed
    monkeypatch.setattr(ascatsc_r.shutil, "which", lambda _: None)


def result_folder(tmp_path: Path, viewer: str | None = ".cnvscope.h5ad", older=False) -> Path:
    d = tmp_path / "ascat_results"
    d.mkdir()
    rda = d / "result_object.Rda"
    rda.write_bytes(b"RDX3\n not a real R object")
    if viewer:
        f = d / f"result_object{viewer}"
        if viewer.endswith(".h5ad"):
            load(FIX).write(f)
        else:
            shutil.copy(FIX, f)
        if older:
            t = rda.stat().st_mtime - 3600
            os.utime(f, (t, t))
    return d


@pytest.mark.parametrize("viewer", [".cnvscope.h5ad", ".h5ad", ".cnvscope.json.gz"])
def test_rda_opens_from_viewer_file(tmp_path, no_r, viewer):
    d = result_folder(tmp_path, viewer)
    ds = load(d / "result_object.Rda")
    assert ds.adata.n_obs == load(FIX).adata.n_obs


def test_older_viewer_file_is_ignored(tmp_path, no_r):
    d = result_folder(tmp_path, older=True)
    assert ascatsc_r.viewer_file(d / "result_object.Rda") is None
    with pytest.raises(RuntimeError, match="needs R"):
        load(d / "result_object.Rda")


def test_without_viewer_file_r_is_used(tmp_path, no_r):
    d = result_folder(tmp_path, viewer=None)
    with pytest.raises(RuntimeError, match="needs R"):
        load(d / "result_object.Rda")


def run(*args):
    return CliRunner().invoke(main, ["--mode", "dna", "--tool", "ascatsc", *args, "--no-view"])


def test_cli_folder_and_input_file_use_viewer_file(tmp_path, no_r):
    d = result_folder(tmp_path)
    # the folder: picks the viewer file itself
    r = run("--input", str(d), "--outdir", str(tmp_path / "o1"))
    assert r.exit_code == 0, r.output
    assert "result_object.cnvscope.h5ad" in r.output
    # the .Rda named: opened from the viewer file beside it
    r = run("--input-file", str(d / "result_object.Rda"), "--outdir", str(tmp_path / "o2"))
    assert r.exit_code == 0, r.output
    assert (tmp_path / "o2" / "result_object.h5ad").exists()


def test_cli_force_reads_the_object_itself(tmp_path, no_r):
    d = result_folder(tmp_path)
    r = run(
        "--input-file", str(d / "result_object.Rda"), "--outdir", str(tmp_path / "o"), "--force"
    )
    assert r.exit_code != 0 and "needs R" in r.output
