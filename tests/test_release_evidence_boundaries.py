"""Release evidence must survive relocation and fail closed at each boundary."""
import copy
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from benchmarks.compare import release_provenance as p
from benchmarks.compare import release_validation as rv


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value))


def seal_attempt(destination, record):
    history = destination / 'attempt_history' / record['cell'] / 'test-attempt.json'
    state = {'cell': record['cell'], 'status': 'success', 'attempt': str(history),
             'report': p.fingerprint(destination / 'results' / (record['cell'] + '.json'))}
    write_json(history, state)
    write_json(destination / 'attempts' / (record['cell'] + '.json'), state)


@pytest.fixture
def legacy(tmp_path, monkeypatch):
    data = tmp_path / 'external data'
    data.mkdir()
    paths = {key: data / name for key, name in (
        ('ref_gff', 'ref.gff3'), ('ref_fa', 'ref.fa'), ('tgt_fa', 'target.fa'))}
    annotation = ('##gff-version 3\n'
                  'chr1\ttest\tgene\t1\t30\t.\t+\t.\tID=g1\n'
                  'chr1\ttest\tmRNA\t1\t30\t.\t+\t.\tID=tx1;Parent=g1\n'
                  'chr1\ttest\texon\t1\t30\t.\t+\t.\tID=e1;Parent=tx1\n'
                  'chr1\ttest\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=tx1\n')
    paths['ref_gff'].write_text(annotation)
    for key in ('ref_fa', 'tgt_fa'):
        paths[key].write_text('>chr1\nATG' + 'AAA' * 8 + 'TAA\n')
    inputs = {name: {**p.fingerprint(paths[key]), 'size_bytes': paths[key].stat().st_size,
                     'fingerprint_status': 'complete', 'changed_during_hash': False}
              for key, name in (('ref_gff', 'reference_annotation'), ('ref_fa', 'reference_genome'),
                                ('tgt_fa', 'target_genome'))}
    source = tmp_path / 'legacy'
    record = {'cell': 'demo', 'arms': {}}
    for arm in (rv.NEW_LABEL, rv.OLD_LABEL):
        output = source / 'cells' / 'demo' / arm / (arm + '.gff3')
        output.parent.mkdir(parents=True)
        output.write_text(annotation)
        write_json(output.parent / 'lifton_output/run_manifest.json',
                   {'inputs': inputs, 'run': {'status': 'success'}})
        table = output.parent.parent / 'eval' / (arm + '.transcripts.tsv')
        table.parent.mkdir(exist_ok=True)
        with table.open('w') as stream:
            writer = csv.DictWriter(stream, delimiter='\t', fieldnames=[
                'ref_mrna_id', 'is_coding', 'recovered', 'protein_identity', 'status'])
            writer.writeheader()
            writer.writerow(dict(ref_mrna_id='tx1', is_coding=1, recovered=1, protein_identity=1, status='ok'))
        record['arms'][arm] = {'output': str(output)}
    write_json(source / 'results/demo.json', record)
    monkeypatch.setattr(rv, 'HERE', tmp_path / 'relocated source' / 'benchmarks' / 'compare')
    monkeypatch.setattr(rv.fc, '_full_paths', lambda bid: (_ for _ in ()).throw(
        AssertionError('worktree-relative input lookup is forbidden for legacy records')))
    monkeypatch.setattr(rv, '_evaluation_evidence', lambda: {'test': 'fixed'})
    return source, tmp_path / 'rescored', paths


def test_legacy_rescore_uses_preserved_inputs_from_relocated_source(legacy):
    source, destination, paths = legacy
    before = {f: f.read_bytes() for f in source.rglob('*') if f.is_file()}
    result = rv.rescore('demo', source, destination, log=lambda _: None)
    assert result['inputs']['ref_gff'] == p.fingerprint(paths['ref_gff'])
    assert result['common_set']['n'] == 1
    assert result['gate']['provenance_verified'] is False
    assert all(f.read_bytes() == body for f, body in before.items())


@pytest.mark.parametrize('damage', ['missing', 'inconsistent', 'modified'])
def test_legacy_rescore_rejects_missing_inconsistent_or_changed_inputs(legacy, damage):
    source, destination, paths = legacy
    manifest = source / 'cells/demo' / rv.OLD_LABEL / 'lifton_output/run_manifest.json'
    document = json.loads(manifest.read_text())
    if damage == 'missing':
        del document['inputs']['reference_annotation']
    elif damage == 'inconsistent':
        document['inputs']['reference_annotation']['sha256'] = 'a' * 64
    else:
        paths['ref_gff'].write_text('changed bytes')
    write_json(manifest, document)
    with pytest.raises((ValueError, RuntimeError), match='[Ii]nput'):
        rv.rescore('demo', source, destination, log=lambda _: None)
    assert not (destination / 'results/demo.json').exists()


def test_run_cell_evaluation_preserves_source_sidecars(legacy, monkeypatch):
    source, destination, paths = legacy
    monkeypatch.setattr(rv.fc, '_full_paths', lambda bid: paths)
    monkeypatch.setattr(rv.fc, '_bench', lambda bid: {'species': 'synthetic', 'cross_species': False})
    sidecars = [Path(str(paths['ref_gff']) + suffix) for suffix in ('_db', '.eval_db')]
    sidecars += [Path(str(paths[key]) + '.fai') for key in ('ref_fa', 'tgt_fa')]
    for path in sidecars:
        path.write_text('chr1\t30\t6\t30\t31\n' if path.suffix == '.fai' else 'preserve previous evidence')
    def arm(bid, label, tree, paths, anndb, config, log):
        historical_label = rv.NEW_LABEL if label == 'candidate' else rv.OLD_LABEL
        output = source / 'cells/demo' / historical_label / (historical_label + '.gff3')
        sidecar = Path(str(output) + '_db')
        sidecar.write_text('preserve previous evidence')
        sidecars.append(sidecar)
        receipt = {'validation': rv.validate_output(output), 'wall_seconds': 1, 'peak_rss_mb': 1,
                   'evidence': {'runtime': {'version': '1.0.13' if label == 'candidate' else '1.0.12'},
                                'source': {'commit': 'a' * 40}}}
        write_json(output.parent / 'completion.json', receipt)
        return output, receipt
    monkeypatch.setattr(rv, '_run_arm', arm)
    result = rv.run_cell('demo', rv.Configuration(destination, source, source, threads=1), log=lambda _: None)
    assert set(result['arms']) == {'candidate', 'reference'}
    assert result['arms']['candidate']['version'] == '1.0.13'
    assert result['arms']['reference']['commit'] == 'a' * 40
    assert all(path.read_text() == ('chr1\t30\t6\t30\t31\n' if path.suffix == '.fai'
                                    else 'preserve previous evidence') for path in sidecars)


def report(output, bid='demo'):
    from tests.test_release_validation import _record
    record = rv._finish(_record())
    record.update(schema_version=2, cell=bid)
    write_json(output / 'results' / (bid + '.json'), record)
    return record


@pytest.mark.parametrize('expected', [None, []])
def test_merge_requires_explicit_nonempty_expected_cells(tmp_path, expected):
    report(tmp_path)
    assert rv.merge(tmp_path, expected) == 1
    assert not json.loads((tmp_path / 'release_validation.json').read_text())['gate_pass']


@pytest.mark.parametrize('damage', ['duplicate', 'unexpected', 'failed_retry', 'stale_output', 'unsealed'])
def test_merge_rejects_ambiguous_or_stale_reports(tmp_path, damage):
    record = report(tmp_path)
    if damage == 'duplicate':
        write_json(tmp_path / 'results/copy.json', record)
    elif damage == 'unexpected':
        report(tmp_path, 'unrequested')
    elif damage == 'failed_retry':
        write_json(tmp_path / 'failures/demo.json', {'cell': 'demo', 'error': 'failed retry'})
    elif damage == 'stale_output':
        artifact = tmp_path / 'output.gff3'
        artifact.write_text('original')
        for details in record['arms'].values():
            details['output_artifact'] = p.fingerprint(artifact)
        artifact.write_text('changed')
        write_json(tmp_path / 'results/demo.json', record)
    assert rv.merge(tmp_path, ['demo']) == 1


def test_receipt_rejects_native_tool_change_during_run(tmp_path, monkeypatch):
    source = {'root': str(tmp_path), 'commit': 'a' * 40}
    monkeypatch.setattr(p, 'snapshot', lambda _: source)
    monkeypatch.setattr(p, 'runtime', lambda *args: {'tools': {'minimap2': 'changed'}})
    output, manifest = tmp_path / 'out.gff3', tmp_path / 'manifest.json'
    output.write_text('##gff-version 3\n')
    write_json(manifest, {'run': {'status': 'success'}})
    expected = {'source': source, 'inputs': {}, 'argv': ['python'], 'environment': {},
                'cwd': str(tmp_path), 'runtime': {'tools': {'minimap2': 'original'}}}
    with pytest.raises(RuntimeError, match='[Rr]untime|[Tt]ool'):
        p.write_receipt(tmp_path / 'receipt.json', expected=expected, output=output,
                        manifest=manifest, profile=SimpleNamespace(exit_code=0, wall_clock_seconds=1, peak_rss_mb=1),
                        validation={'complete': True})


def test_evaluator_evidence_detects_same_version_dependency_changes(tmp_path, monkeypatch):
    package = tmp_path / 'package.py'
    package.write_text('original = 1\n')
    dist = SimpleNamespace(version='1.0', files=[Path('package.py')], requires=[], locate_file=lambda name: tmp_path / name)
    import importlib.metadata
    monkeypatch.setattr(importlib.metadata, 'distribution', lambda name: dist)
    before = rv._evaluation_evidence()
    package.write_text('modified = 2\n')
    assert rv._evaluation_evidence() != before


def test_evaluator_evidence_covers_lifton_scoring_sources():
    evidence = rv._evaluation_evidence()
    assert 'lifton/align.py' in evidence['sources']
    assert 'lifton/extract_sequence.py' in evidence['sources']
    assert 'lifton/gff3_validator.py' in evidence['sources']


def test_validator_bounds_nonerrors_without_hiding_late_error(tmp_path):
    output = tmp_path / 'many-warnings.gff3'
    lines = ['##gff-version 3']
    for n in range(100):
        lines += [f'chr1\tLiftOn\tgene\t1\t30\t.\t+\t.\tID=g{n}',
                  f'chr1\tLiftOn\tmRNA\t1\t30\t.\t+\t.\tID=t{n};Parent=g{n}']
    lines[-1] += ';protein_identity=2.0'
    output.write_text('\n'.join(lines) + '\n')
    result = rv.validate_output(output)
    assert result['n_warnings'] >= 100
    assert len([i for i in result['issues'] if i['severity'] != 'ERROR']) <= 20
    assert any(i['feature_id'] == 't99' and i['check'] == 'lifton_protein_identity'
               for i in result['issues'])
    assert rv.new_validity_issues(result, {'issues': []})


def test_reference_validation_does_not_require_lifton_attributes(legacy):
    result = rv.validate_output(legacy[2]['ref_gff'], reference=True)
    assert not any(i['check'].startswith('lifton_') for i in result['issues'])


@pytest.mark.parametrize('recovered', ['0', '1'])
def test_reference_extraction_failure_remains_coding_and_unresolved(legacy, recovered):
    from tests.test_release_validation import _record, row
    source, destination, paths = legacy
    paths['ref_fa'].write_text('>other_chromosome\nATGAAATAA\n')
    ref, _ = rv.evaluator.build_reference(str(paths['ref_gff']), str(paths['ref_fa']), log=lambda _: None)
    assert ref['tx1']['is_coding'] is True
    assert ref['tx1']['prot'] == ''
    inventory = rv.reference_inventory(ref)
    assert inventory['coding_ids'] == ['tx1']
    assert inventory['unresolved_coding_ids'] == ['tx1']
    record = _record()
    rows = [dict(row('tx1', '', recovered=recovered), n_cds_lifted='0', lifted_prot_len='0')]
    record['common_set'] = rv.common_set(rows, rows, expected_ids=ref,
                                        expected_coding_ids=inventory['coding_ids'])
    record['reference_inventory'] = inventory
    assert not rv._finish(record)['gate_pass']
    assert record['gate']['reference_scoring_resolved'] is False


def test_stop_codon_alone_is_not_reference_coding_evidence(legacy):
    _, _, paths = legacy
    paths['ref_gff'].write_text(paths['ref_gff'].read_text().replace('\tCDS\t', '\tstop_codon\t'))
    ref, _ = rv.evaluator.build_reference(str(paths['ref_gff']), str(paths['ref_fa']), log=lambda _: None)
    assert ref['tx1']['is_coding'] is False and ref['tx1']['n_cds'] == 0


@pytest.fixture
def complete_report(legacy, monkeypatch):
    source, destination, paths = legacy
    import sys
    write_json(destination / 'campaign.json', {
        'expected_cells': ['demo'], 'candidate': str(Path.cwd()), 'reference': str(destination),
        'python': sys.executable, 'threads': 8, 'copies': False, 'rescore': None})
    record = rv.rescore('demo', source, destination, log=lambda _: None)
    record.pop('legacy_source')
    record['reference_inventory'] = {'schema_version': 1, 'all_ids': ['tx1'], 'coding_ids': ['tx1'],
                                     'unresolved_coding_ids': []}
    record['protocol'] = {'fresh_alignment': True, 'copies': False}
    source_evidence = {'root': str(source), 'commit': 'a' * 40}
    monkeypatch.setattr(p, 'snapshot', lambda _: source_evidence)
    runtime = lambda python, root, cwd, env: {
        'tools': 'fixed', 'version': '1.0.13' if rv.NEW_LABEL in cwd else '1.0.12'}
    monkeypatch.setattr(p, 'runtime', runtime)
    for arm, details in record['arms'].items():
        output = Path(details['output'])
        expected = {'source': source_evidence, 'inputs': record['inputs'],
                    'argv': ['python'], 'cwd': str(output.parent), 'environment': {},
                    'runtime': runtime('python', source, str(output.parent), {})}
        write_json(output.parent / 'expected.json', expected)
        receipt_path = output.parent / 'completion.json'
        p.write_receipt(receipt_path, expected=expected, output=output,
                        manifest=output.parent / 'lifton_output/run_manifest.json',
                        profile=SimpleNamespace(exit_code=0, wall_clock_seconds=1, peak_rss_mb=1),
                        validation=details['validity'])
        details['receipt'] = p.fingerprint(receipt_path)
        details['provenance_verified'] = True
    rv._finish(record)
    rv._write_report(destination, 'demo', record)
    seal_attempt(destination, record)
    assert rv.merge(destination, ['demo']) == 0
    return destination, record


@pytest.mark.parametrize('damage', ['campaign', 'attempt', 'history', 'history_changed'])
def test_required_campaign_and_attempt_evidence_cannot_disappear(complete_report, damage):
    destination, _ = complete_report
    paths = {'campaign': destination / 'campaign.json', 'attempt': destination / 'attempts/demo.json',
             'history': destination / 'attempt_history/demo/test-attempt.json'}
    if damage == 'history_changed':
        state = json.loads(paths['history'].read_text())
        state['status'] = 'running'
        write_json(paths['history'], state)
    else:
        paths[damage].unlink()
    assert rv.merge(destination, ['demo']) == 1
    assert not json.loads((destination / 'release_validation.json').read_text())['gate_pass']


def test_resume_cannot_recreate_deleted_campaign(complete_report):
    destination, _ = complete_report
    (destination / 'campaign.json').unlink()
    with pytest.raises(SystemExit) as exc:
        rv.main(['demo', '--resume', '--run-id', destination.name,
                 '--output-root', str(destination.parent), '--reference-tree', str(destination)])
    assert exc.value.code != 0
    assert not (destination / 'campaign.json').exists()


@pytest.mark.parametrize('damage', ['output', 'table', 'receipt', 'report', 'evaluator',
                                    'failed_retry', 'running_retry', 'duplicate', 'unexpected',
                                    'duplicate_expected', 'missing_expected', 'empty_expected'])
def test_verified_report_cannot_mask_later_invalid_evidence(complete_report, monkeypatch, damage):
    destination, record = complete_report
    expected = ['demo']
    if damage in ('output', 'table', 'receipt'):
        key = {'output': 'output_artifact', 'table': 'transcript_table', 'receipt': 'receipt'}[damage]
        Path(record['arms'][rv.NEW_LABEL][key]['path']).write_text('changed')
    elif damage == 'report':
        path = destination / 'results/demo.json'
        path.write_text(path.read_text() + '\n')
    elif damage == 'evaluator':
        monkeypatch.setattr(rv, '_evaluation_evidence', lambda: {'test': 'changed'})
    elif damage in ('failed_retry', 'running_retry'):
        write_json(destination / 'attempts/demo.json', {'status': damage.split('_')[0]})
    elif damage == 'duplicate':
        write_json(destination / 'results/copy.json', record)
    elif damage == 'unexpected':
        other = copy.deepcopy(record)
        other['cell'] = 'other'
        rv._write_report(destination, 'other', other)
    elif damage == 'duplicate_expected':
        expected += ['demo']
    elif damage == 'missing_expected':
        expected = None
    else:
        expected = []
    assert rv.merge(destination, expected) == 1
    assert not json.loads((destination / 'release_validation.json').read_text())['gate_pass']


def test_cli_failed_retry_invalidates_previous_success(complete_report, monkeypatch):
    destination, record = complete_report
    def fail(*args):
        raise RuntimeError('retry failed')
    monkeypatch.setattr(rv, 'run_cell', fail)
    assert rv.main(['demo', '--resume', '--run-id', destination.name,
                    '--output-root', str(destination.parent), '--reference-tree', str(destination),
                    '--candidate-tree', str(Path.cwd())]) == 1
    assert rv.merge(destination, ['demo']) == 1
    assert len(list((destination / 'attempt_history/demo').glob('*.json'))) == 2
    assert json.loads((destination / 'results/demo.json').read_text())['gate_pass'] is True


def test_cli_rejects_changed_campaign_expected_set(complete_report):
    destination, record = complete_report
    write_json(destination / 'campaign.json', {'expected_cells': ['demo', 'missing']})
    assert rv.merge(destination, ['demo']) == 1


def test_successful_retry_supersedes_failure_without_erasing_history(complete_report, monkeypatch):
    destination, record = complete_report
    failed = {'cell': 'demo', 'status': 'failed', 'error': 'previous failure'}
    history = destination / 'attempt_history/demo/previous-failed.json'
    write_json(history, failed)
    write_json(destination / 'failures/demo.json', failed)
    write_json(destination / 'attempts/demo.json', failed)
    monkeypatch.setattr(rv, 'run_cell', lambda *args: rv._write_report(destination, 'demo', record))
    assert rv.main(['demo', '--resume', '--run-id', destination.name,
                    '--output-root', str(destination.parent), '--reference-tree', str(destination),
                    '--candidate-tree', str(Path.cwd())]) == 0
    assert json.loads(history.read_text()) == failed
    assert rv.merge(destination, ['demo']) == 0


def test_coding_membership_cannot_be_reclassified_by_the_evaluator():
    from tests.test_release_validation import _record, row
    record = _record()
    noncoding = dict(row('a'), is_coding='0')
    record['common_set'] = rv.common_set([noncoding], [noncoding], expected_ids=['a'], expected_coding_ids=['a'])
    assert not rv._finish(record)['gate_pass']
    assert record['gate']['coding_membership_consistent'] is False


def test_corrupt_report_cannot_leave_stale_successful_summary(complete_report):
    destination, _ = complete_report
    (destination / 'results/demo.json').write_text('broken JSON')
    with pytest.raises(ValueError):
        rv.merge(destination, ['demo'])
    assert not json.loads((destination / 'release_validation.json').read_text())['gate_pass']


def test_validator_error_count_without_all_error_identities_is_unresolved():
    from tests.test_release_validation import _record
    record = _record()
    record['arms'][rv.NEW_LABEL]['validity']['n_errors'] = 1
    assert not rv._finish(record)['gate_pass']


def test_evaluator_change_during_reference_validation_refuses_report(legacy, monkeypatch):
    source, destination, _ = legacy
    state = {'source': 'original'}
    monkeypatch.setattr(rv, '_evaluation_evidence', lambda: dict(state))
    validate = rv.validate_output
    def changing_validator(path, *, reference=False):
        result = validate(path, reference=reference)
        if reference:
            state['source'] = 'changed'
        return result
    monkeypatch.setattr(rv, 'validate_output', changing_validator)
    with pytest.raises(RuntimeError, match='Evaluator'):
        rv.rescore('demo', source, destination, log=lambda _: None)
    assert not (destination / 'results/demo.json').exists()


@pytest.mark.parametrize('changed', ['source', 'runtime'])
def test_merge_rechecks_live_arm_provenance(complete_report, monkeypatch, changed):
    destination, _ = complete_report
    if changed == 'source':
        monkeypatch.setattr(p, 'snapshot', lambda _: {'root': 'changed', 'commit': 'b' * 40})
    else:
        monkeypatch.setattr(p, 'runtime', lambda *args: {'tools': 'changed'})
    assert rv.merge(destination, ['demo']) == 1


def test_version_neutral_reports_preserve_release_metadata(complete_report):
    destination, record = complete_report
    record['schema_version'] = 3
    record['arms'] = {role: {**record['arms'][label], 'role': role,
                             'version': '1.0.13' if role == 'candidate' else '1.0.12',
                             'commit': 'a' * 40}
                      for role, label in [('candidate', rv.NEW_LABEL), ('reference', rv.OLD_LABEL)]}
    rv._write_report(destination, 'demo', record)
    seal_attempt(destination, record)
    assert rv.merge(destination, ['demo']) == 0
    merged = json.loads((destination / 'release_validation.json').read_text())['records'][0]
    assert set(merged['arms']) == {'candidate', 'reference'}
    assert merged['arms']['candidate']['version'] == '1.0.13'


def test_empty_evaluation_cannot_pass():
    from tests.test_release_validation import _record
    record = _record()
    record['common_set'] = rv.common_set([], [])
    assert not rv._finish(record)['gate_pass']


def test_duplicate_unrecovered_reference_rows_cannot_pass():
    from tests.test_release_validation import _record, row
    record = _record()
    record['common_set'] = rv.common_set([row('a'), row('b', recovered='0'), row('b', recovered='0')],
                                         [row('a'), row('b', recovered='0')])
    assert not rv._finish(record)['gate_pass']


def test_missing_reference_evaluation_rows_cannot_pass():
    from tests.test_release_validation import _record, row
    record = _record()
    record['common_set'] = rv.common_set([row('a')], [row('a'), row('b', recovered='0')])
    assert not rv._finish(record)['gate_pass']


def test_replaced_report_is_archived_byte_for_byte(tmp_path):
    record = report(tmp_path)
    path = tmp_path / 'results/demo.json'
    original = path.read_bytes()
    digest = p.sha256(path)
    rv._write_report(tmp_path, 'demo', record)
    assert (tmp_path / 'report_history/demo' / (digest + '.json')).read_bytes() == original


def test_summary_separates_presence_cds_scores_and_unresolved():
    from tests.test_release_validation import row
    rows = [dict(row('scored'), n_cds_lifted='1', lifted_prot_len='10'),
            dict(row('unresolved', ''), n_cds_lifted='1', lifted_prot_len='10'),
            dict(row('no_cds', ''), n_cds_lifted='0', lifted_prot_len='0')]
    result = rv._summarize(rows, {})
    assert result['n_present_coding'] == 3
    assert result['n_cds_recovered_coding'] == 2
    assert result['n_scored_coding'] == 1
    assert result['n_unresolved_coding'] == 1


def test_legacy_table_change_while_reading_refuses_report(legacy, monkeypatch):
    source, destination, _ = legacy
    original = rv.gene_level.read_transcript_tsv
    def changing_reader(path):
        rows = original(path)
        Path(path).write_text('changed while reading')
        return rows
    monkeypatch.setattr(rv.gene_level, 'read_transcript_tsv', changing_reader)
    with pytest.raises((ValueError, RuntimeError), match='[Tt]able|[Aa]rtifact'):
        rv.rescore('demo', source, destination, log=lambda _: None)
    assert not (destination / 'results/demo.json').exists()


def test_role_version_must_match_recorded_runtime(complete_report):
    destination, record = complete_report
    record['schema_version'] = 3
    record['arms'] = {role: {**record['arms'][label], 'role': role,
                             'version': 'invented', 'commit': 'a' * 40}
                      for role, label in [('candidate', rv.NEW_LABEL), ('reference', rv.OLD_LABEL)]}
    rv._write_report(destination, 'demo', record)
    assert rv.merge(destination, ['demo']) == 1


def test_dependency_evidence_follows_runtime_closure_and_requested_extras(tmp_path, monkeypatch):
    import importlib.metadata
    package = tmp_path / 'package.py'
    package.write_text('value = 1\n')
    def distribution(name):
        requires = {'gffutils': ['argh[cli]>=1', 'pytest; extra == "test"',
                                'wrong-platform; python_version < "2"'],
                    'argh': ['argcomplete; extra == "cli"']}.get(name, [])
        return SimpleNamespace(version='1.0', files=[Path('package.py')], requires=requires,
                               locate_file=lambda name: tmp_path / name)
    monkeypatch.setattr(importlib.metadata, 'distribution', distribution)
    evidence = p.dependency_evidence()
    assert 'argh' in evidence and 'argcomplete' in evidence
    assert 'pytest' not in evidence and 'wrong-platform' not in evidence


def test_an_absent_optional_extra_is_recorded_not_fatal(tmp_path, monkeypatch):
    """`mappy` is an optional extra and CI asserts it is ABSENT.

    Listing it beside the runtime requirements made every release-evidence test
    fail on a correctly provisioned machine -- a supported configuration read as
    a broken one. Its absence is recorded so the evidence still says which of
    the two environments produced it.
    """
    import importlib.metadata
    package = tmp_path / 'package.py'
    package.write_text('value = 1\n')

    def distribution(name):
        if name == 'mappy':
            raise importlib.metadata.PackageNotFoundError(name)
        return SimpleNamespace(version='1.0', files=[Path('package.py')], requires=[],
                               locate_file=lambda name: tmp_path / name)

    monkeypatch.setattr(importlib.metadata, 'distribution', distribution)
    evidence = p.dependency_evidence()
    assert evidence['mappy'] == {'installed': False}
    assert evidence['numpy']['version'] == '1.0'


def test_an_installed_optional_extra_is_fingerprinted(tmp_path, monkeypatch):
    import importlib.metadata
    package = tmp_path / 'package.py'
    package.write_text('value = 1\n')
    monkeypatch.setattr(
        importlib.metadata, 'distribution',
        lambda name: SimpleNamespace(version='2.28', files=[Path('package.py')],
                                     requires=[],
                                     locate_file=lambda name: tmp_path / name))
    evidence = p.dependency_evidence()
    assert evidence['mappy']['version'] == '2.28'
    assert 'sha256' in evidence['mappy']


def test_missing_mandatory_dependency_fails_evidence(tmp_path, monkeypatch):
    import importlib.metadata
    def missing(name):
        raise importlib.metadata.PackageNotFoundError(name)
    monkeypatch.setattr(importlib.metadata, 'distribution', missing)
    with pytest.raises(RuntimeError, match='Missing required dependency'):
        p.dependency_evidence()


def test_fresh_cli_requires_explicit_reference_tree(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(rv, 'run_cell', lambda *args: None)
    monkeypatch.setattr(rv, 'merge', lambda *args: 1)
    with pytest.raises(SystemExit) as exc:
        rv.main(['demo', '--run-id', 'new', '--output-root', str(tmp_path)])
    assert exc.value.code == 2
    assert '--reference-tree' in capsys.readouterr().err


def test_identically_truncated_tables_fail_against_reference_inventory():
    from tests.test_release_validation import _record, row
    record = _record()
    record['common_set'] = rv.common_set([row('a')], [row('a')], expected_ids={'a', 'b'})
    assert record['common_set']['missing_evaluation_new_ids'] == ['b']
    assert record['common_set']['missing_evaluation_old_ids'] == ['b']
    assert not rv._finish(record)['gate_pass']


def test_evaluation_cannot_substitute_an_unexpected_reference_id():
    from tests.test_release_validation import _record, row
    record = _record()
    record['common_set'] = rv.common_set([row('a'), row('extra')], [row('a'), row('extra')],
                                         expected_ids={'a', 'b'})
    assert record['common_set']['unexpected_evaluation_new_ids'] == ['extra']
    assert not rv._finish(record)['gate_pass']


def _legacy_distribution(tmp_path, listings):
    """A distribution whose file list is its SOURCES.txt -- `setup.py`, not
    installed, beside `package.py`, installed. `listings` are the other
    inventory files it has."""
    (tmp_path / 'package.py').write_text('value = 1\n')
    texts = {'SOURCES.txt': 'setup.py\npackage.py\n'}
    texts.update({name: 'listed\n' for name in listings})
    return SimpleNamespace(
        version='0.2.7', requires=[], files=[Path('setup.py'), Path('package.py')],
        locate_file=lambda name: tmp_path / name, read_text=texts.get)


@pytest.mark.parametrize('listings', [set(), {'installed-files.txt'}])
def test_a_legacy_source_inventory_hashes_what_is_installed(tmp_path, monkeypatch, listings):
    """An egg-info install whose file list is its SOURCES.txt -- the source
    tree. Python 3.10 returns that list even when installed-files.txt exists.
    A legacy interlap install listed nine files never installed and failed
    evidence on an intact env."""
    import importlib.metadata
    monkeypatch.setattr(importlib.metadata, 'distribution',
                        lambda name: _legacy_distribution(tmp_path, listings))
    evidence = p.dependency_evidence()
    assert evidence['interlap']['files'] == 1
    assert evidence['interlap']['inventory'] == 'egg-info SOURCES.txt'


def test_a_file_missing_from_an_install_record_is_fatal(tmp_path, monkeypatch):
    import importlib.metadata
    monkeypatch.setattr(importlib.metadata, 'distribution',
                        lambda name: _legacy_distribution(tmp_path, {'RECORD'}))
    with pytest.raises(RuntimeError, match='lists setup.py as installed'):
        p.dependency_evidence()
