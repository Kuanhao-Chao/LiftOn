"""The release gate must propagate failed stages and preserve prior evidence."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
STUB = r'''
import os, pathlib, shutil, sys
p = pathlib.Path(sys.argv[0])
args = sys.argv[1:]
stage = os.environ.get('FAIL_STAGE', '')
wheel = 'venv' in p.parts
if p.name == 'lifton':
    if '-o' in args:
        out = pathlib.Path(args[args.index('-o') + 1])
        if stage == ('wheel' if wheel else 'lift'):
            sys.exit(17)
        out.parent.mkdir(parents=True, exist_ok=True)
        body = '##gff-version 3\n'
        if stage == 'mismatch' and wheel:
            body += '# different wheel output\n'
        if stage == 'thread_mismatch' and args[args.index('-t') + 1] == '8':
            body += '# different threaded output\n'
        out.write_text(body)
    else:
        print('v1.0.12 --no-orf-stop-completion --intermediate-dir')
elif p.name == 'gff3-validate':
    if stage == 'validate' and '-h' not in args:
        sys.exit(19)
    print('VALID')
elif p.name == 'pip':
    assert '--no-cache-dir' in args
    if stage == 'install':
        sys.exit(21)
elif '-m' in args and args[args.index('-m') + 1] == 'build':
    if stage == 'build':
        sys.exit(23)
    out = pathlib.Path(args[args.index('--outdir') + 1])
    out.mkdir(parents=True, exist_ok=True)
    (out / 'lifton-1.0.12-py3-none-any.whl').write_text('wheel')
    (out / 'lifton-1.0.12.tar.gz').write_text('sdist')
elif '-m' in args and args[args.index('-m') + 1] == 'venv':
    out = pathlib.Path(args[-1]) / 'bin'
    out.mkdir(parents=True, exist_ok=True)
    for name in ('python', 'pip', 'lifton', 'gff3-validate'):
        shutil.copy2(p, out / name)
'''


@pytest.fixture
def release_workspace(tmp_path):
    root = tmp_path / 'source with spaces'
    (root / 'benchmarks').mkdir(parents=True)
    (root / 'test').mkdir()
    script = root / 'benchmarks' / 'release_gates.sh'
    shutil.copyfile(ROOT / 'benchmarks' / 'release_gates.sh', script)
    bin_dir = tmp_path / 'tools'
    bin_dir.mkdir()
    for name in ('python', 'pip', 'lifton', 'gff3-validate'):
        path = bin_dir / name
        path.write_text(f'#!{sys.executable}\n' + STUB)
        path.chmod(0o755)
    return script, bin_dir / 'python', tmp_path / 'evidence'


def _run(work, stage=''):
    script, python, out = work
    env = {**os.environ, 'LIFTON_PY': str(python), 'FAIL_STAGE': stage}
    return subprocess.run(['bash', str(script), str(out)], env=env,
                          text=True, capture_output=True, timeout=30)


@pytest.mark.parametrize('stage', [
    'lift', 'validate', 'build', 'install', 'wheel', 'mismatch', 'thread_mismatch',
])
def test_failed_stage_is_a_failed_gate(release_workspace, stage):
    result = _run(release_workspace, stage)
    assert result.returncode != 0, result.stdout + result.stderr
    assert '[DONE] exit=0' not in result.stdout


def test_success_requires_all_stages(release_workspace):
    result = _run(release_workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    assert '[DONE] exit=0' in result.stdout
    assert (release_workspace[2] / 'wheel_chr22' / 'lifton.gff3').is_file()


def test_existing_evidence_is_never_deleted(release_workspace):
    out = release_workspace[2]
    out.mkdir()
    evidence = out / 'previous-result.txt'
    evidence.write_text('keep me')
    result = _run(release_workspace)
    assert result.returncode != 0
    assert evidence.read_text() == 'keep me'
