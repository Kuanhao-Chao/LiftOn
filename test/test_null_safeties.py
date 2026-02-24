import pytest
from lifton import lifton_utils
from lifton import run_miniprot

def test_miniprot_id_mapping_null_safety():
    # Pass None explicitly to simulate miniprot crash producing an empty file
    ref_dict, m_dict = lifton_utils.miniprot_id_mapping(None)
    
    # Assert it returns empty dictionaries instead of crashing
    assert isinstance(ref_dict, dict)
    assert len(ref_dict) == 0
    
    assert isinstance(m_dict, dict)
    assert len(m_dict) == 0

def test_process_miniprot_null_safety():
    # Attempting to process a miniprot transcript when feature_db is None
    # Provide dummy mocks for other args
    lifton_gene = run_miniprot.process_miniprot(
        mtrans=None,
        ref_db=None,
        m_feature_db=None,
        tree_dict={},
        tgt_fai=None,
        ref_proteins=None,
        ref_trans=None,
        ref_features_dict=None,
        fw_score=None,
        m_id_2_ref_id_trans_dict={},
        ref_features_len_dict={},
        ref_trans_exon_num_dict={},
        ref_features_reverse_dict={},
        args=None
    )
    
    # process_miniprot should explicitly return early giving None
    assert lifton_gene is None
