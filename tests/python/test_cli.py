"""The one-command run: cnvscope --mode --tool --input --outdir."""

import json
from pathlib import Path

import anndata as ad
import pandas as pd
from click.testing import CliRunner
from test_tool_readers import write_chisel

from cnvscope.cli import main

FIX = Path(__file__).parents[1] / "fixtures" / "ascatsc_json"


def run(*args):
    return CliRunner().invoke(main, list(args), catch_exceptions=False)


def test_chisel_run_writes_outputs(tmp_path):
    write_chisel(tmp_path / "chisel_out")
    out = tmp_path / "results"
    r = run(
        "--mode",
        "dna",
        "--tool",
        "chisel",
        "--input",
        str(tmp_path / "chisel_out"),
        "--outdir",
        str(out),
        "--no-view",
    )
    assert r.exit_code == 0, r.output
    a = ad.read_h5ad(out / "chisel_out.h5ad")
    assert a.n_obs == 3 and a.uns["tool"] == "CHISEL"
    cells = pd.read_csv(out / "chisel_out_cells.tsv", sep="\t")
    assert cells["cell"].tolist() == ["AAAC", "GGTT", "TTAG"] and "ploidy" in cells
    meta = json.loads((out / "chisel_out_run.json").read_text())
    assert meta["mode"] == "dna" and meta["tool"] == "chisel" and meta["cells"] == 3


def test_shorthand_flag_and_ascatsc_folder(tmp_path):
    r = run(
        "--mode",
        "dna",
        "--ascatsc",
        "--input",
        str(FIX),
        "--outdir",
        str(tmp_path),
        "--name",
        "demo",
        "--no-view",
    )
    assert r.exit_code == 0, r.output
    assert (tmp_path / "demo.h5ad").exists()
    assert "demo_export.json.gz" in r.output  # found inside the folder


def test_helpful_errors(tmp_path):
    d = str(tmp_path)
    r = CliRunner().invoke(main, ["--mode", "rna", "--tool", "chisel", "--input", d, "--outdir", d])
    assert r.exit_code != 0 and "CHISEL is a scDNA tool" in r.output and "numbat" in r.output
    r = CliRunner().invoke(main, ["--mode", "rna", "--tool", "numbat", "--input", d, "--outdir", d])
    assert "not supported yet" in r.output
    r = CliRunner().invoke(main, ["--mode", "dna", "--tool", "chisel", "--input", d])
    assert "missing --outdir" in r.output
    r = CliRunner().invoke(main, ["--mode", "dna", "--tool", "nope", "--input", d, "--outdir", d])
    assert "unknown tool" in r.output
    r = CliRunner().invoke(
        main, ["--mode", "dna", "--chisel", "--hiscanner", "--input", d, "--outdir", d]
    )
    assert "choose one tool" in r.output
    r = CliRunner().invoke(main, ["--mode", "dna", "--tool", "chisel", "--input", d, "--outdir", d])
    assert "does not look like CHISEL output" in r.output
    r = CliRunner().invoke(
        main, ["--mode", "multiome", "--tool", "chisel", "--input", d, "--outdir", d]
    )
    assert "multiome" in r.output


def test_help_and_tools():
    assert "cnvscope --mode dna --tool chisel" in run().output
    assert "chisel" in run("tools").output and "planned" in run("tools").output


def test_busy_port_moves_to_the_next_free_one():
    import socket

    from cnvscope.cli import free_port

    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        s.listen()
        busy = s.getsockname()[1]
        assert free_port("127.0.0.1", busy) > busy


def test_input_file_inside_the_folder_or_as_a_path(tmp_path):
    # a name inside --input: the dataset is named after the folder
    r = run(
        "--mode",
        "dna",
        "--ascatsc",
        "--input",
        str(FIX),
        "--input-file",
        "demo_export.json.gz",
        "--outdir",
        str(tmp_path),
        "--no-view",
    )
    assert r.exit_code == 0, r.output
    assert (tmp_path / "ascatsc_json.h5ad").exists()
    assert "also found" not in r.output
    # a full path on its own (other spelling): named after the file
    r = run(
        "--mode",
        "dna",
        "--tool",
        "ascat.sc",
        "--input_file",
        str(FIX / "demo_export.json.gz"),
        "--outdir",
        str(tmp_path),
        "--no-view",
    )
    assert r.exit_code == 0, r.output
    assert (tmp_path / "demo_export.h5ad").exists()


def test_input_file_not_found_lists_the_folder(tmp_path):
    r = CliRunner().invoke(
        main,
        [
            "--mode",
            "dna",
            "--ascatsc",
            "--input",
            str(FIX),
            "--input-file",
            "results.rds",
            "--outdir",
            str(tmp_path),
        ],
    )
    assert r.exit_code != 0
    assert "results.rds not found" in r.output and "demo_export.json.gz" in r.output


def test_second_run_reuses_the_result_until_the_input_changes(tmp_path):
    src = tmp_path / "chisel_out"
    write_chisel(src)
    out = tmp_path / "results"
    args = ("--mode", "dna", "--chisel", "--input", str(src), "--outdir", str(out), "--no-view")
    assert "reading CHISEL output" in run(*args).output
    second = run(*args).output
    assert "unchanged since the last run" in second and "reading CHISEL" not in second
    assert "reading CHISEL output" in run(*args, "--force").output
    # change the input: it is read again
    calls = src / "calls" / "calls.tsv"
    calls.write_text(calls.read_text() + "")
    import os

    os.utime(calls, ns=(calls.stat().st_atime_ns, calls.stat().st_mtime_ns + 10**9))
    assert "reading CHISEL output" in run(*args).output
