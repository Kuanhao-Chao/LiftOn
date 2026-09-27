"""Qualification uses the existing controller and immutable per-attempt evidence."""
import copy
import json
from pathlib import Path
import sys

import pytest

from benchmarks.compare import build_controller as controller
from benchmarks.compare import release_validation as rv


@pytest.fixture
def config(tmp_path):
    inputs = {}
    for name in ('ref_gff', 'ref_fa', 'tgt_fa'):
        path = tmp_path / name
        path.write_text('input\n')
        inputs[name] = str(path)
    return {'schema_version': 1, 'python': sys.executable,
            'candidate': {'root': str(tmp_path / 'candidate'), 'sha': 'a' * 40},
            'reference': {'root': str(tmp_path / 'reference'), 'sha': 'b' * 40},
            'cells': {'demo': {'inputs': inputs, 'species': 'synthetic', 'cross_species': False,
                               'annotation_database': 'RefSeq', 'copies': False,
                               'candidate_mode': 'safe', 'reference_mode': 'safe', 'full_job': True}}}


def test_controller_dispatches_qualification_with_resource_accounting(tmp_path, config):
    from benchmarks.compare import release_qualification as q
    frozen = {'configuration': config, 'sources': {}, 'inputs': {'demo': {}}}
    cells = controller.build_cells('qualification', ['demo'], run_dir=tmp_path, policy=controller.Policy(),
                                   dataset_registry=tmp_path / 'datasets.json', qualification=frozen)
    cell, = cells
    assert cell['kind'] == 'release_qualification'
    assert cell['command'] == [sys.executable, '-m', 'benchmarks.compare.release_qualification',
                               '--cell', str(Path(cell['cell_dir']) / 'cell.json')]
    assert cell['threads'] == 8 and cell['scheduler_thread_cost'] == 16
    assert cell['full_job'] is True
    resources = {'available_gib': 300, 'load1': 0}
    assert controller.launch_allowed(cell, [cell], resources, controller.Policy())[0]
    assert not controller.launch_allowed(cell, [cell, cell], resources, controller.Policy())[0]
    assert controller._hard_timeout_seconds(cell, controller.Policy()) == 2 * 8 * 3600
    controller.validate_plan_layout({'run_dir': str(tmp_path), 'cells': cells})
    assert q.attempt_output(cell, 1) != q.attempt_output(cell, 2)


@pytest.mark.parametrize('damage', ['duplicate', 'missing', 'unknown', 'unsafe'])
def test_qualification_selection_requires_exact_unambiguous_ids(config, damage):
    from benchmarks.compare import release_qualification as q
    requested = {'duplicate': ['demo', 'demo'], 'missing': [], 'unknown': ['absent'], 'unsafe': ['../demo']}[damage]
    with pytest.raises(ValueError):
        q.select_ids(config, requested)


def test_configuration_freezes_absolute_inputs_and_exact_sources(tmp_path, config, monkeypatch):
    from benchmarks.compare import release_qualification as q
    from benchmarks.compare import release_provenance as p
    def snapshot(root):
        role = Path(root).name
        return {'root': str(root), 'commit': config[role]['sha']}
    monkeypatch.setattr(p, 'snapshot', snapshot)
    monkeypatch.setattr(p, 'runtime', lambda *args: {'fixed': True})
    monkeypatch.setattr(p, 'evaluation_evidence', lambda: {'fixed': True})
    frozen = q.prepare(config, ['demo'])
    assert frozen['inputs']['demo']['ref_gff']['sha256']
    Path(config['cells']['demo']['inputs']['ref_gff']).write_text('changed\n')
    assert q.prepare(config, ['demo']) != frozen
    bad = copy.deepcopy(config)
    bad['candidate']['sha'] = 'not-exact'
    with pytest.raises(ValueError):
        q.prepare(bad, ['demo'])


def test_qualification_uses_explicit_inputs_without_worktree_lookup(tmp_path, config, monkeypatch):
    from benchmarks.compare import release_qualification as q
    monkeypatch.setattr(rv.fc, '_full_paths', lambda _: pytest.fail('must use pinned paths'))
    frozen = {'configuration': config, 'sources': {}, 'inputs': {'demo': {}}}
    cell = q.build_cell('demo', tmp_path / 'cell', 8, frozen)
    run_config = q.run_configuration(cell, 1)
    assert run_config.paths == {key: Path(path) for key, path in config['cells']['demo']['inputs'].items()}
    assert run_config.metadata['species'] == 'synthetic'
    assert run_config.output == q.attempt_output(cell, 1)
    assert run_config.copies is False


def test_qualification_retry_keeps_failed_directory(tmp_path, config, monkeypatch):
    from benchmarks.compare import release_qualification as q
    frozen = {'configuration': config, 'sources': {}, 'inputs': {'demo': {}}}
    cell = q.build_cell('demo', tmp_path / 'cell', 8, frozen)
    monkeypatch.setattr(q, 'verify_frozen', lambda _: None)
    def failure(*args, **kwargs):
        raise RuntimeError('native failure')
    monkeypatch.setattr(rv, 'run_cell', failure)
    assert q.run_attempt(cell, 1) == 1
    failed = q.attempt_output(cell, 1) / 'attempts/demo.json'
    original = failed.read_bytes()
    assert json.loads(original)['status'] == 'failed'
    assert q.run_attempt(cell, 2) == 1
    assert failed.read_bytes() == original
    assert (q.attempt_output(cell, 2) / 'attempts/demo.json').exists()
    with pytest.raises(FileExistsError):
        q.run_attempt(cell, 1)


def test_controller_qualification_plan_roundtrip_and_tampering(tmp_path, config, monkeypatch):
    from benchmarks.compare import release_qualification as q
    from benchmarks.compare import release_provenance as p
    monkeypatch.setattr(p, 'snapshot', lambda root: {'root': root, 'commit': config[Path(root).name]['sha']})
    monkeypatch.setattr(p, 'runtime', lambda *args: {'fixed': True})
    monkeypatch.setattr(p, 'evaluation_evidence', lambda: {'fixed': True})
    base_provenance = {'test': 'controller'}
    base_provenance['fingerprint'] = controller.canonical_hash(base_provenance)
    monkeypatch.setattr(controller, 'collect_provenance', lambda **kwargs: dict(base_provenance))
    run_dir, plan = controller.create_plan(
        run_id='test', stage='qualification', requested_ids=['demo'], runs_root=tmp_path / 'runs',
        repo_root=tmp_path, registry=tmp_path / 'benchmarks.json', dataset_registry=tmp_path / 'datasets.json',
        baseline=tmp_path / 'baseline.json', policy=controller.Policy(), qualification=config)
    controller.initialize_run(run_dir, plan)
    assert controller.load_plan(run_dir) == plan
    controller.assert_matching_provenance(plan)
    changed = copy.deepcopy(plan)
    changed['cells'][0]['command'].append('--unexpected')
    changed['cells'][0]['fingerprint'] = controller._cell_fingerprint(changed['cells'][0], changed['provenance']['fingerprint'])
    changed['fingerprint'] = controller._plan_fingerprint(changed)
    with pytest.raises(ValueError, match='execution'):
        controller.validate_plan_integrity(changed)
    Path(config['cells']['demo']['inputs']['tgt_fa']).write_text('changed target\n')
    with pytest.raises(RuntimeError, match='provenance'):
        controller.assert_matching_provenance(plan)


@pytest.mark.parametrize('policy', [controller.Policy(threads_per_cell=16), controller.Policy(min_available_gib=128),
                                    controller.Policy(max_worker_threads=64), controller.Policy(max_full=3)])
def test_qualification_cannot_exceed_resource_envelope(tmp_path, config, policy):
    with pytest.raises(ValueError, match='approved'):
        controller.build_cells('qualification', ['demo'], run_dir=tmp_path, policy=policy,
                               dataset_registry=tmp_path / 'datasets.json', qualification={'configuration': config})
