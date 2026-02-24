import pytest
import copy
from lifton import run_liftoff
from lifton.lifton_class import Lifton_GENE, LiftOn_FEATURE

class DummyFeatureDB:
    def __init__(self, children_dict):
        self.children_dict = children_dict
        
    def children(self, locus, featuretype=None, level=1, order_by=None):
        if locus.id in self.children_dict:
            res = self.children_dict[locus.id]
            if featuretype:
                res = [f for f in res if f.featuretype == featuretype]
            return res
        return []

class DummyLocus:
    def __init__(self, id, featuretype):
        self.id = id
        self.featuretype = featuretype
        self.attributes = {"ID": [id]}
        self.source = "LiftOn"
        self.start = 100
        self.end = 200
        self.seqid = "chr1"

class DummyArgs:
    def __init__(self):
        self.annotation_database = "GENCODE"
        self.no_orf_search = True
        self.debug = False

def test_process_liftoff_recursion_preserves_root_gene():
    # Build a simulated structure where a transcript lacks exons natively
    # Gene -> ncRNA -> ncRNA_nested (simulating an intermediate structure)
    gene_locus = DummyLocus("gene1", "gene")
    tx_locus = DummyLocus("tx1", "transcript")
    intermediate_locus = DummyLocus("tx1_nested", "ncRNA")
    
    # Feature DB logic
    children_dict = {
        "gene1": [tx_locus],
        "tx1": [intermediate_locus],
        "tx1_nested": [] # No exons beneath this
    }
    mock_db = DummyFeatureDB(children_dict)
    
    # Provide a mocked lifton_gene that represents the root
    # Note: we are directly calling the recursion segment, so we bypass ENTRY_FEATURE=True init
    root_gene = Lifton_GENE("ref_gene1", gene_locus, {}, {}, {}, DummyArgs())
    root_gene.entry.id = "gene1"
    
    ref_db = DummyFeatureDB({})
    ref_id_2_m_id_trans_dict = {}
    m_feature_db = DummyFeatureDB({})
    tree_dict = {}
    ref_proteins = {}
    ref_trans = {}
    ref_features_dict = {}
    fw_score = None
    fw_chain = None
    args = DummyArgs()
    tgt_fai = {}

    # The buggy process_liftoff used to return 'lifton_gene = process_liftoff(...)' for child loops,
    # overwriting the root_gene. Now we confirm it correctly traverses down AND returns the same root object unmodified.
    
    returned_gene = run_liftoff.process_liftoff(
        lifton_gene=root_gene, 
        locus=gene_locus, 
        ref_db=ref_db, 
        l_feature_db=mock_db, 
        ref_id_2_m_id_trans_dict=ref_id_2_m_id_trans_dict, 
        m_feature_db=m_feature_db, 
        tree_dict=tree_dict, 
        tgt_fai=tgt_fai, 
        ref_proteins=ref_proteins, 
        ref_trans=ref_trans, 
        ref_features_dict=ref_features_dict, 
        fw_score=fw_score, 
        fw_chain=fw_chain, 
        args=args, 
        ENTRY_FEATURE=False
    )
    
    # It should still be the ROOT GENE, not a LiftOn_FEATURE
    assert isinstance(returned_gene, Lifton_GENE)
    assert returned_gene.ref_gene_id == "ref_gene1"
