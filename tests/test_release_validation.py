"""Release comparisons must not hide losses, missing scores, or missing runs."""
import copy
import math

import pytest

from benchmarks.compare import release_validation as rv


def row(identifier, pi='0.9', *, recovered='1', status='ok'):
    return {'ref_mrna_id': identifier, 'protein_identity': pi,
            'is_coding': '1', 'recovered': recovered, 'status': status}


def test_net_gain_does_not_hide_lost_ids():
    result = rv.common_set([row('b'), row('c')], [row('a')])
    assert result['lost_ids'] == ['a']
    assert result['added_ids'] == ['b', 'c']
    assert result['n'] == 0


@pytest.mark.parametrize('missing', ['', 'NaN', 'inf', '-inf', None])
def test_missing_or_nonfinite_scores_are_explicit(missing):
    result = rv.common_set([row('a', missing)], [row('a')])
    assert result['lost_ids'] == []
    assert result['newly_unscored_ids'] == ['a']
    assert result['n'] == 0
    assert result['mean_new'] is None


def test_shared_unscored_is_not_claimed_as_compared():
    result = rv.common_set([row('a', '', status='map_failed')],
                           [row('a', '', status='map_failed')])
    assert result['shared_recovered'] == 1
    assert result['n'] == 0
    assert result['unscored_new_ids'] == ['a']
    assert result['statuses_new'] == {'map_failed': 1}


def test_duplicate_reference_rows_are_not_overwritten_silently():
    result = rv.common_set([row('a'), row('a', '1.0')], [row('a')])
    assert result['duplicate_new_ids'] == ['a']


def test_low_identity_additions_do_not_dilute_common_set():
    result = rv.common_set([row('a'), row('b', '0.4')], [row('a')])
    assert result['n'] == 1
    assert result['n_regressed'] == 0
    assert result['delta'] == 0
    assert math.isfinite(result['mean_new'])


def test_validity_errors_cannot_cancel_each_other():
    old = {'issues': [{'check': 'containment', 'feature_id': 'a', 'severity': 'ERROR'}]}
    new = {'issues': [{'check': 'containment', 'feature_id': 'b', 'severity': 'ERROR'}]}
    assert rv.new_validity_issues(new, old) == [['containment', 'b', 1]]


def _record():
    arm = {'recall': {'n_recovered_coding': 1, 'mean_protein_identity': 0.9},
           'completed': True, 'provenance_verified': True,
           'validity': {'exit': 0, 'n_errors': 0, 'n_warnings': 0,
                        'valid': True, 'issues': [], 'complete': True}}
    return {'arms': {rv.NEW_LABEL: copy.deepcopy(arm), rv.OLD_LABEL: copy.deepcopy(arm)},
            'common_set': rv.common_set([row('a')], [row('a')])}


def test_known_complete_comparison_passes():
    assert rv._finish(_record())['gate_pass']


@pytest.mark.parametrize('key,value', [
    ('completed', False), ('provenance_verified', False), ('validity', {}),
])
def test_missing_evidence_cannot_pass(key, value):
    record = _record()
    record['arms'][rv.NEW_LABEL][key] = value
    assert not rv._finish(record)['gate_pass']


def test_unscored_models_require_an_explicit_resolution():
    record = _record()
    record['common_set'] = rv.common_set([row('a', '')], [row('a', '')])
    assert not rv._finish(record)['gate_pass']


def test_missing_cds_is_resolved_but_not_a_scored_protein():
    missing = {**row('a', '', status='map_failed'), 'n_cds_lifted': '0', 'lifted_prot_len': '0'}
    record = _record()
    record['common_set'] = rv.common_set([missing], [missing])
    assert record['common_set']['n'] == 0
    assert record['common_set']['missing_cds_new_ids'] == ['a']
    assert rv._finish(record)['gate_pass']


def test_new_missing_cds_is_a_coding_model_regression():
    missing = {**row('a', '', status='map_failed'), 'n_cds_lifted': '0', 'lifted_prot_len': '0'}
    record = _record()
    record['common_set'] = rv.common_set([missing], [row('a')])
    assert record['common_set']['lost_coding_model_ids'] == ['a']
    assert not rv._finish(record)['gate_pass']


def test_cli_no_cells_is_not_success(tmp_path):
    with pytest.raises(SystemExit) as exc:
        rv.main([])
    assert exc.value.code != 0
