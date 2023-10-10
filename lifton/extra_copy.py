from lifton import lifton_class
import copy
from intervaltree import Interval, IntervalTree

def find_extra_copy(m_feature_db, tree_dict, aa_id_2_m_id_dict, gene_info_dict, trans_info_dict, trans_2_gene_dict, gene_copy_num_dict, trans_copy_num_dict, fw):

    EXTRA_COPY_MINIPROT_COUNT = 0
    NEW_LOCUS_MINIPROT_COUNT = 0

    for mtrans in m_feature_db.features_of_type('mRNA'):
        chromosome = mtrans.seqid
        mtrans_id = mtrans.attributes["ID"][0]
        mtrans_interval = Interval(mtrans.start, mtrans.end, mtrans_id)
        ovps = tree_dict[chromosome].overlap(mtrans_interval)

        # Add new transcript loci into interval tree
        gene_interval = Interval(mtrans.start, mtrans.end, mtrans.attributes["ID"])
        tree_dict[chromosome].add(gene_interval)

        if len(ovps) == 0:

            extra_cp_trans_id = aa_id_2_m_id_dict[mtrans_id]
            gene_entry_base = copy.deepcopy(mtrans)
            trans_entry_base = copy.deepcopy(mtrans)
            # Find the extra copy for know gene
            if extra_cp_trans_id in trans_info_dict.keys():
                EXTRA_COPY_MINIPROT_COUNT += 1
                extra_cp_gene_id = trans_2_gene_dict[extra_cp_trans_id] 

                #######################################
                # Step 5.1: Create the gene entry
                #######################################
                gene_attrs = gene_info_dict[extra_cp_gene_id]
                # print("gene_attrs : ", gene_attrs)
                Lifton_gene_ecp = lifton_class.Lifton_GENE(gene_entry_base)
                new_extra_cp_gene_id = Lifton_gene_ecp.update_gene_info(extra_cp_gene_id, chromosome, mtrans.start, mtrans.end, gene_attrs, gene_copy_num_dict)


                #######################################
                # Step 5.2: Create the transcript entry
                #######################################
                trans_attrs = trans_info_dict[extra_cp_trans_id]
                # print("trans_attrs: ", trans_attrs)
                new_extra_cp_trans_id = Lifton_gene_ecp.create_new_transcript(False, extra_cp_gene_id, extra_cp_trans_id, trans_entry_base, chromosome, mtrans.start, mtrans.end, trans_attrs, gene_copy_num_dict, trans_copy_num_dict)

                #######################################
                # Step 5.3: Create the exon entry
                #######################################

                #######################################
                # Step 5.4: Create the CDS entry
                #######################################
                cdss = m_feature_db.children(mtrans, featuretype='CDS')  # Replace 'exon' with the desired child feature type
                # print("cdss len: ", len(cdss))
                for cds in list(cdss):
                    # entry.attributes["ID"]
                    Lifton_gene_ecp.add_exon(new_extra_cp_trans_id, cds)
                    cds_copy = copy.deepcopy(cds)
                    Lifton_gene_ecp.add_cds(new_extra_cp_trans_id, cds_copy)
                Lifton_gene_ecp.write_entry(fw)
            else:
                # Reason it's missing => the mRNA does not belong to gene (vdj segments)
                NEW_LOCUS_MINIPROT_COUNT += 1
                #######################################
                # Step 5.1: Create the gene entry
                #######################################
                extra_cp_gene_id = "gene-LiftOn"
                Lifton_gene_ecp = lifton_class.Lifton_GENE(gene_entry_base)
                gene_attrs = lifton_class.Lifton_GENE_info({}, extra_cp_gene_id)
                new_extra_cp_gene_id = Lifton_gene_ecp.update_gene_info(extra_cp_gene_id, chromosome, mtrans.start, mtrans.end, gene_attrs, gene_copy_num_dict)

                #######################################
                # Step 5.2: Create the transcript entry
                #######################################
                trans_attrs = lifton_class.Lifton_TRANS_info({}, extra_cp_trans_id, extra_cp_gene_id)
                new_extra_cp_trans_id = Lifton_gene_ecp.create_new_transcript(True, extra_cp_gene_id, extra_cp_trans_id, trans_entry_base, chromosome, mtrans.start, mtrans.end, trans_attrs, gene_copy_num_dict, trans_copy_num_dict)

                #######################################
                # Step 5.3: Create the exon entry
                #######################################
                #######################################
                # Step 5.4: Create the CDS entry
                #######################################
                cdss = m_feature_db.children(mtrans, featuretype='CDS')  # Replace 'exon' with the desired child feature type
                # print("cdss len: ", len(cdss))
                for cds in list(cdss):
                    Lifton_gene_ecp.add_exon(extra_cp_trans_id, cds)
                    cds_copy = copy.deepcopy(cds)
                    Lifton_gene_ecp.add_cds(extra_cp_trans_id, cds_copy)
                Lifton_gene_ecp.write_entry(fw)

    return EXTRA_COPY_MINIPROT_COUNT, NEW_LOCUS_MINIPROT_COUNT
