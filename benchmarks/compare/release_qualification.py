"""Release qualification adapter for build_controller; no separate scheduler."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys

from . import release_provenance as provenance
from . import release_validation as validation
from lifton.run_manifest import atomic_write_json


MODES = {'safe': (), 'stream': ('--stream',), 'inmemory': ('--inmemory-liftoff',),
         'stream-inmemory': ('--stream', '--inmemory-liftoff')}


def select_ids(configuration, requested=None):
    available = configuration.get('cells', {})
    ids = list(available) if requested is None else list(requested)
    if (not ids or len(ids) != len(set(ids)) or not isinstance(available, dict)
            or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', value) or value not in available for value in ids)):
        raise ValueError('Qualification requires a nonempty unique set of known, path-safe cells')
    return ids


def prepare(configuration, ids):
    """Normalize and content-address the selected source, runtime and inputs."""
    ids = select_ids(configuration, ids)
    configuration = json.loads(json.dumps(configuration))
    if configuration.get('schema_version') != 1:
        raise ValueError('Unsupported qualification configuration schema')
    python = Path(configuration['python'])
    if not python.is_absolute() or not python.is_file():
        raise ValueError('Qualification python must be an absolute interpreter path')
    sources, runtimes, inputs = {}, {}, {}
    for role in ('candidate', 'reference'):
        source = configuration[role]
        if not re.fullmatch('[0-9a-f]{40}', source['sha']) or not Path(source['root']).is_absolute():
            raise ValueError(f'{role} requires an absolute source path and exact lowercase commit SHA')
        sources[role] = provenance.snapshot(source['root'])
        if sources[role]['commit'] != source['sha']:
            raise ValueError(f'{role} source does not match the selected commit')
        source['root'] = sources[role]['root']
        runtimes[role] = provenance.runtime(str(python), source['root'], source['root'],
                                            provenance.isolated_env(source['root']))
    for identifier in ids:
        case = configuration['cells'][identifier]
        if set(case['inputs']) != {'ref_gff', 'ref_fa', 'tgt_fa'}:
            raise ValueError(f'{identifier}: all three explicit inputs are required')
        for name in ('copies', 'full_job', 'cross_species'):
            if not isinstance(case.get(name), bool):
                raise ValueError(f'{identifier}: {name} must be an explicit boolean')
        for role in ('candidate', 'reference'):
            if case.get(role + '_mode') not in MODES:
                raise ValueError(f'{identifier}: unsupported {role} execution mode')
        if case.get('annotation_database') not in ('RefSeq', 'Ensembl', 'GENCODE', 'CHESS'):
            raise ValueError(f'{identifier}: unsupported annotation_database')
        if not isinstance(case.get('species'), str) or not case['species']:
            raise ValueError(f'{identifier}: species metadata is required')
        if any(not Path(path).is_absolute() for path in case['inputs'].values()):
            raise ValueError(f'{identifier}: input paths must be absolute')
        inputs[identifier] = {name: provenance.fingerprint(path) for name, path in case['inputs'].items()}
        case['inputs'] = {name: record['path'] for name, record in inputs[identifier].items()}
    configuration['cells'] = {identifier: configuration['cells'][identifier] for identifier in ids}
    return {'configuration': configuration, 'expected_cells': ids, 'sources': sources,
            'runtimes': runtimes, 'inputs': inputs, 'evaluator': provenance.evaluation_evidence()}


def verify_frozen(frozen):
    if prepare(frozen['configuration'], frozen['expected_cells']) != frozen:
        raise RuntimeError('Qualification source, inputs, runtime or evaluator changed')


def validate_policy(policy):
    if (policy.threads_per_cell > 8 or policy.max_full > 2 or policy.max_worker_threads > 32
            or policy.min_available_gib < 256):
        raise ValueError('Qualification exceeds the approved 8-thread/cell, 2-full, 32-thread, 256-GiB-reserve policy')
    cost = max(2 * policy.threads_per_cell, policy.scheduler_threads_per_cell or 0)
    if cost > policy.max_worker_threads:
        raise ValueError('Qualification cell plus native thread reserve exceeds the thread budget')
    return cost


def build_cell(identifier, cell_dir, threads, frozen):
    case = frozen['configuration']['cells'][identifier]
    cell_dir = Path(cell_dir)
    return {'id': 'qualification__' + identifier, 'kind': 'release_qualification', 'benchmark': identifier,
            'mode': 'qualification', 'threads': threads, 'full_job': case['full_job'],
            'command': [sys.executable, '-m', 'benchmarks.compare.release_qualification',
                        '--cell', str(cell_dir / 'cell.json')],
            'environment': {'PYTHONPATH': str(Path(__file__).resolve().parents[2])},
            'qualification': frozen, 'cell_dir': str(cell_dir),
            'artifacts': {'result_json': str(cell_dir / 'qualification_result.json')}}


def attempt_output(cell, attempt):
    if not isinstance(attempt, int) or isinstance(attempt, bool) or attempt < 1:
        raise ValueError('Qualification attempt must be a positive integer')
    return Path(cell['cell_dir']) / 'qualification_attempts' / f'{attempt:04d}'


def run_configuration(cell, attempt):
    config = cell['qualification']['configuration']
    case = config['cells'][cell['benchmark']]
    return validation.Configuration(
        attempt_output(cell, attempt), Path(config['candidate']['root']), Path(config['reference']['root']),
        python=config['python'], threads=cell['threads'], copies=case['copies'],
        paths={name: Path(path) for name, path in case['inputs'].items()}, metadata=case,
        role_options={role: MODES[case[role + '_mode']] for role in ('candidate', 'reference')})


def run_attempt(cell, attempt):
    verify_frozen(cell['qualification'])
    config = run_configuration(cell, attempt)
    config.output.mkdir(parents=True, exist_ok=False)
    bid = cell['benchmark']
    atomic_write_json(config.output / 'campaign.json', {
        'expected_cells': [bid], 'controller_cell': cell['id'], 'controller_attempt': attempt,
        'controller_fingerprint': cell.get('fingerprint'), 'qualification': cell['qualification']})
    history = config.output / 'attempt_history' / bid / 'execution.json'
    active = config.output / 'attempts' / (bid + '.json')
    state = {'cell': bid, 'status': 'running', 'attempt': str(history)}
    atomic_write_json(history, state)
    atomic_write_json(active, state)
    try:
        validation.run_cell(bid, config)
        verify_frozen(cell['qualification'])
        state.update(status='success', report=provenance.fingerprint(config.output / 'results' / (bid + '.json')))
    except Exception as error:
        state.update(status='failed', error=str(error))
        atomic_write_json(config.output / 'failures' / (bid + '.json'), state)
    atomic_write_json(history, state)
    atomic_write_json(active, state)
    result = validation.merge(config.output, [bid])
    atomic_write_json(cell['artifacts']['result_json'], {
        'schema_version': 1, 'cell': cell['id'], 'attempt': attempt, 'output': str(config.output),
        'fingerprint': cell.get('fingerprint'),
        'summary': provenance.fingerprint(config.output / 'release_validation.json'),
        'gate_pass': state['status'] == 'success' and result == 0})
    return 0 if state['status'] == 'success' and result == 0 else 1


def validate_artifacts(cell, started_ns):
    """Read-only verification used by controller completion, resume and audit."""
    from . import build_controller as controller
    errors, artifacts = [], {}
    try:
        result_path = Path(cell['artifacts']['result_json'])
        if not controller._artifact_is_fresh(result_path, started_ns):
            raise ValueError('Qualification result is missing, empty or stale')
        result = json.loads(result_path.read_text())
        status = json.loads((Path(cell['cell_dir']) / 'status.json').read_text())
        attempt = status['attempts']
        output = attempt_output(cell, attempt)
        if (result.get('schema_version') != 1 or result.get('cell') != cell['id']
                or result.get('attempt') != attempt or result.get('fingerprint') != cell['fingerprint']
                or result.get('output') != str(output) or result.get('gate_pass') is not True):
            raise ValueError('Qualification result does not match the current controller attempt')
        if result.get('summary') != provenance.fingerprint(output / 'release_validation.json'):
            raise ValueError('Qualification summary changed')
        summary = json.loads((output / 'release_validation.json').read_text())
        bid = cell['benchmark']
        if (summary.get('gate_pass') is not True or summary.get('expected_cells') != [bid]
                or len(summary.get('records', [])) != 1):
            raise ValueError('Qualification summary is incomplete')
        report_path = output / 'results' / (bid + '.json')
        report = json.loads(report_path.read_text())
        if report != summary['records'][0]:
            raise ValueError('Qualification summary and sealed report disagree')
        campaign = json.loads((output / 'campaign.json').read_text())
        if (campaign.get('qualification') != cell['qualification'] or campaign.get('controller_attempt') != attempt
                or campaign.get('controller_cell') != cell['id']
                or campaign.get('controller_fingerprint') != cell['fingerprint']):
            raise ValueError('Qualification campaign differs from immutable controller configuration')
        if not validation._finish(report)['gate_pass']:
            raise ValueError('Qualification gates no longer pass')
        config = run_configuration(cell, attempt)
        if (report.get('threads') != config.threads or report['protocol'].get('copies') != config.copies
                or report['protocol'].get('role_options') != {k: list(v) for k, v in config.role_options.items()}):
            raise ValueError('Qualification report execution modes differ from the controller configuration')
        errors.extend(validation._report_errors(output, report_path, report, provenance.evaluation_evidence()))
        # Seal only evidence, not logs, DB indexes or temporary evaluator inputs.
        evidence = [result_path, output / 'release_validation.json', output / 'campaign.json', report_path,
                    output / 'report_receipts' / (bid + '.json'), output / 'attempts' / (bid + '.json'),
                    output / 'attempt_history' / bid / 'execution.json']
        for role, arm in report['arms'].items():
            receipt = Path(arm['receipt']['path'])
            receipt_data = json.loads(receipt.read_text())
            if receipt_data['evidence']['source'] != cell['qualification']['sources'][role]:
                raise ValueError(f'{role} executed a different source from the controller configuration')
            expected_argv = validation._argv(
                config.python, config.paths, config.metadata['annotation_database'], config.threads,
                Path(arm['output']), root=cell['qualification']['configuration'][role]['root'],
                copies=config.copies, options=config.role_options[role])
            if receipt_data['evidence']['argv'] != expected_argv:
                raise ValueError(f'{role} executed a different command from the controller configuration')
            evidence.extend((Path(arm['output']), Path(arm['transcript_table']['path']), receipt,
                             receipt.parent / 'expected.json', Path(receipt_data['manifest']['path'])))
        artifacts = {str(path.relative_to(Path(cell['cell_dir']))): controller._success_artifact_record(path)
                     for path in evidence}
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as error:
        errors.append(str(error))
    return errors, {'artifacts': artifacts}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cell', type=Path, required=True)
    args = parser.parse_args(argv)
    from . import build_controller as controller
    run_dir = args.cell.resolve().parents[2]
    plan = controller.load_plan(run_dir)
    cell = controller._cell_for(plan, args.cell.parent.name)
    if args.cell != Path(cell['cell_dir']) / 'cell.json' or cell['kind'] != 'release_qualification':
        raise ValueError('Cell does not match immutable controller plan')
    status = controller._read_status(cell)
    if status.get('state') != 'running':
        raise ValueError('Qualification must be launched by the running controller cell')
    return run_attempt(cell, status['attempts'])


if __name__ == '__main__':
    sys.exit(main())
