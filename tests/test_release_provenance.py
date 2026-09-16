"""A cached filename or old successful process cannot stand in for evidence."""
import json
from types import SimpleNamespace

import pytest

from benchmarks.compare import release_provenance as p


@pytest.fixture
def receipt(tmp_path, monkeypatch):
    source = {'root': str(tmp_path), 'commit': 'a' * 40}
    monkeypatch.setattr(p, 'snapshot', lambda root: source)
    output = tmp_path / 'out.gff3'
    output.write_text('##gff-version 3\n')
    manifest = tmp_path / 'manifest.json'
    manifest.write_text(json.dumps({'run': {'status': 'success'}}))
    expected = {'source': source, 'inputs': {}, 'argv': ['lifton'], 'runtime': {'version': '1'},
                'cwd': str(tmp_path), 'environment': {}}
    monkeypatch.setattr(p, 'runtime', lambda *args: {'version': '1'})
    path = tmp_path / 'completion.json'
    profile = SimpleNamespace(exit_code=0, wall_clock_seconds=1, peak_rss_mb=2)
    p.write_receipt(path, expected=expected, output=output, manifest=manifest,
                    profile=profile, validation={'complete': True})
    return path, expected, output, manifest


def test_complete_matching_receipt_can_resume(receipt):
    assert p.read_receipt(*receipt)['wall_seconds'] == 1


@pytest.mark.parametrize('changed', ['output', 'manifest', 'command', 'tool', 'source'])
def test_changed_artifact_or_provenance_cannot_resume(receipt, changed):
    path, expected, output, manifest = receipt
    if changed == 'output':
        output.write_text('different bytes')
    elif changed == 'manifest':
        manifest.write_text(json.dumps({'run': {'status': 'failed'}}))
    elif changed == 'command':
        expected['argv'] += ['--stream']
    elif changed == 'tool':
        expected['runtime']['version'] = '2'
    else:
        expected['source']['commit'] = 'b' * 40
    assert p.read_receipt(path, expected, output, manifest) is None


def test_nonempty_output_without_receipt_cannot_resume(receipt):
    path, expected, output, manifest = receipt
    path.unlink()
    assert output.stat().st_size > 0
    assert p.read_receipt(path, expected, output, manifest) is None


def test_partial_manifest_cannot_receive_receipt(receipt):
    path, expected, output, manifest = receipt
    manifest.write_text(json.dumps({'run': {'status': 'partial'}}))
    with pytest.raises(RuntimeError, match='Incomplete/partial'):
        p.write_receipt(path, expected=expected, output=output, manifest=manifest,
                        profile=SimpleNamespace(exit_code=0), validation={'complete': True})


@pytest.mark.parametrize('in_tree', [False, True])
def test_actual_cli_bootstrap_rejects_editable_submodule_fallthrough(tmp_path, in_tree):
    import os
    from pathlib import Path
    import subprocess
    import sys
    from benchmarks.compare import release_validation as rv
    root, external = tmp_path / 'candidate', tmp_path / 'external'
    package = root / 'lifton'
    package.mkdir(parents=True)
    external.mkdir()
    (package / '__init__.py').write_text('__version__ = "test"\n')
    (package / 'lifton.py').write_text('def main():\n    from lifton import late\n    print(late.VALUE)\n')
    (external / 'late.py').write_text('VALUE = "wrong tree"\n')
    (external / 'sitecustomize.py').write_text(
        'import importlib.util, sys\n'
        'class Fallback:\n'
        '    @classmethod\n'
        '    def find_spec(cls, fullname, path=None, target=None):\n'
        '        if fullname == "lifton.late":\n'
        f'            return importlib.util.spec_from_file_location(fullname, {str(external / "late.py")!r})\n'
        'sys.meta_path.append(Fallback)\n')
    if in_tree:
        (package / 'late.py').write_text('VALUE = "correct tree"\n')
    argv = rv._argv(sys.executable, dict(ref_gff='a', ref_fa='b', tgt_fa='c'),
                    'RefSeq', 1, tmp_path / 'out.gff3', root=root)
    result = subprocess.run(argv, cwd=tmp_path,
                            env=dict(os.environ, PYTHONPATH=os.pathsep.join((str(root), str(external)))),
                            capture_output=True, text=True, timeout=30)
    if in_tree:
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == 'correct tree'
    else:
        assert result.returncode != 0
        assert 'late' in result.stderr and 'wrong tree' not in result.stdout
        control = list(argv)
        control[2] = 'from lifton.lifton import main; main()'
        unguarded = subprocess.run(control, cwd=tmp_path,
                                  env=dict(os.environ, PYTHONPATH=os.pathsep.join((str(root), str(external)))),
                                  capture_output=True, text=True, timeout=30)
        assert unguarded.returncode == 0 and unguarded.stdout.strip() == 'wrong tree'


def test_guard_rejects_an_already_loaded_outside_module(tmp_path):
    import subprocess
    import sys
    source = ('import sys, types\n'
              'm = types.ModuleType("lifton.escape")\n'
              'm.__file__ = "/outside/escape.py"\n'
              'sys.modules["lifton.escape"] = m\n')
    code = source + p.guarded_code(tmp_path, 'print("accepted")')
    result = subprocess.run([sys.executable, '-c', code], cwd=tmp_path,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode != 0 and 'outside' in result.stderr
