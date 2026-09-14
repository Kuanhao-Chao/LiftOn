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
    expected = {'source': source, 'inputs': {}, 'argv': ['lifton'], 'runtime': {'version': '1'}}
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
