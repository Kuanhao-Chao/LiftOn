import shutil
import subprocess

import pytest


def test_miniprot_short_input_smoke(tmp_path):
    executable = shutil.which("miniprot")
    if executable is None:
        pytest.skip("miniprot is not installed")

    target = tmp_path / "target.fa"
    proteins = tmp_path / "proteins.fa"
    target.write_text(">target\nACGTACGTACGT\n")
    proteins.write_text(">prot\nMMMAAA***\n")

    proc = subprocess.run(
        [executable, "--gff-only", str(target), str(proteins)],
        capture_output=True,
        text=True,
        check=False,
    )

    assert proc.returncode == 0
    assert proc.stdout.startswith("##gff-version 3\n")
