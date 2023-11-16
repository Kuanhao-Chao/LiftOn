from lifton import lifton_utils, lifton_class, align, get_id_fraction, variants, extra_copy
import copy, os
from Bio.Seq import Seq
from intervaltree import Interval, IntervalTree

class Lifton_ORF:
    def __init__(self, start, end):
        self.start = start
        self.end = end

class Lifton_Status:
    def __init__(self):
        self.liftoff = 0
        self.miniprot = 0
        self.lifton = 0
        self.annotation = None
        self.status = []

class Lifton_Alignment:
    def __init__(self, extracted_identity, cds_children, alignment_query, alignment_comp, alignment_ref, cdss_protein_boundary, cdss_protein_aln_boundary, extracted_seq, reference_seq, db_entry):
        self.identity = extracted_identity
        self.cds_children = cds_children
        self.query_aln = alignment_query
        self.comp = alignment_comp
        self.ref_aln = alignment_ref
        self.cdss_protein_boundaries = cdss_protein_boundary
        self.cdss_protein_aln_boundaries = cdss_protein_aln_boundary
        self.query_seq = extracted_seq
        self.ref_seq = reference_seq
        self.db_entry = db_entry

    def write_alignment(self, outdir, tool_name, mutation, trans_id):
        outdir_tool = outdir+"/"+tool_name+"/"+mutation+"/"
        os.makedirs(outdir_tool, exist_ok=True)
        outfile= outdir_tool+trans_id+".fa"
        with open(outfile, "w") as fw:
            fw.write("> Reference\n")
            fw.write(self.ref_aln + "\n")
            fw.write("> Target\n")
            fw.write(self.query_aln + "\n")

    
class Lifton_feature:
    def __init__(self, id):
        self.id = id
        self.copy_num = 0
        self.children = set()

class Lifton_GENE:
    def __init__(self, ref_gene_id, gffutil_entry_gene, ref_gene_attrs, tree_dict, ref_features_dict, holder = False):
        ###########################
        # Assigning the reference gene & attributes
        ###########################
        self.entry = gffutil_entry_gene
        self.entry.source = "LiftOn"
        self.entry.featuretype = "gene"
        self.transcripts = {}
        self.copy_num = self.entry.attributes["extra_copy_number"][0] if 'extra_copy_number' in self.entry.attributes else "0"
        # This is used for retrieving the copy num & attributes info
        self.ref_gene_id = ref_gene_id

        if holder:
            # ref_gene_id_in, ref_trans_id_in = lifton_utils.get_ref_ids_liftoff(ref_features_dict, self.entry.id, None)
            attributes = {}
            # if "Parent" in self.entry.attributes.keys():
            attributes["ID"] = self.entry.attributes["Parent"] if "Parent" in self.entry.attributes.keys() else ["LiftOn-gene_" + str(ref_features_dict["LiftOn-gene"].copy_num)]
            ref_features_dict["LiftOn-gene"].copy_num += 1
            self.entry.attributes = attributes
            # print("Holder self.entry.attributes: ", self.entry.attributes)
        else:
            self.entry.attributes = ref_gene_attrs
            self.entry.attributes["ID"] = self.ref_gene_id + "_" + self.copy_num if int(self.copy_num) > 0 else self.ref_gene_id
            if int(self.copy_num) > 0:
                self.entry.attributes["extra_copy_number"] = [self.copy_num]

            # ###########################
            # # Update gene copy number dictionary
            # ###########################
            self.__update_gene_copy(ref_features_dict)
            # print("Nonholder self.entry.attributes: ", self.entry.attributes)
            
        self.entry.id = self.entry.attributes["ID"][0]
        ###########################
        # Adding LiftOn intervals
        ###########################
        gene_interval = Interval(self.entry.start, self.entry.end, self.entry.id)
        if self.entry.seqid not in tree_dict.keys():
            tree_dict[self.entry.seqid] = IntervalTree()
        tree_dict[self.entry.seqid].add(gene_interval)
        

    def __update_gene_copy(self, ref_features_dict):
        if self.ref_gene_id in ref_features_dict.keys():
            ref_features_dict[self.ref_gene_id].copy_num += 1

    def update_gene_info(self, chromosome, start, end):
        self.entry.seqid = chromosome
        self.entry.featuretype = "gene"
        self.entry.start = start
        self.entry.end = end

    def add_miniprot_transcript(self, ref_trans_id, gffutil_entry_trans, ref_trans_attrs, ref_features_dict):
        self.__update_gene_copy(ref_features_dict)
        Lifton_trans = Lifton_TRANS(ref_trans_id, self.ref_gene_id, self.entry.id, self.copy_num, gffutil_entry_trans, ref_trans_attrs)
        self.transcripts[Lifton_trans.entry.id] = Lifton_trans
        return Lifton_trans
    
    def update_trans_info(self, trans_id, chromosome, start, end):
        self.transcripts[trans_id].entry.seqid = chromosome
        self.transcripts[trans_id].entry.featuretype = "mRNA"
        self.transcripts[trans_id].entry.start = start
        self.transcripts[trans_id].entry.end = end

        # transcript_id = lifton_gene.add_transcript(locus, ref_db[ref_trans_id].attributes)
    def add_transcript(self, ref_trans_id, gffutil_entry_trans, ref_trans_attrs):
        Lifton_trans = Lifton_TRANS(ref_trans_id, self.ref_gene_id, self.entry.id, self.copy_num, gffutil_entry_trans, ref_trans_attrs)
        self.transcripts[Lifton_trans.entry.id] = Lifton_trans
        return Lifton_trans

    def remove_transcript(self, transcript_id):
        del self.transcripts[transcript_id]

    def add_exon(self, trans_id, gffutil_entry_exon):
        self.transcripts[trans_id].add_exon(gffutil_entry_exon)

    def add_cds(self, trans_id, gffutil_entry_cds):
        self.transcripts[trans_id].add_cds(gffutil_entry_cds)
                            
    def fix_truncated_protein(self, trans_id, ref_trans_id, fai, ref_proteins, fai_trans, lifton_status):
        ref_protein_seq = ref_proteins[ref_trans_id]
        ref_trans_seq = fai_trans[ref_trans_id]
        if trans_id not in self.transcripts.keys():
            return None, False
        lifton_aln, good_trans = self.transcripts[trans_id].fix_truncated_protein(fai, ref_protein_seq, ref_trans_seq, lifton_status)
        return lifton_aln, good_trans
               
    def update_cds_list(self, trans_id, cds_list):
        self.transcripts[trans_id].update_cds_list(cds_list)
        self.update_boundaries()

    def add_lifton_status_attrs(self, trans_id, lifton_status):
        self.transcripts[trans_id].add_lifton_status_attrs(lifton_status)

    def write_entry(self, fw):
        fw.write(str(self.entry) + "\n")
        for key, trans in self.transcripts.items():
            trans.write_entry(fw)

    def update_boundaries(self):        
        for key, trans in self.transcripts.items():
            self.entry.start = trans.entry.start if trans.entry.start < self.entry.start else self.entry.start
            self.entry.end = trans.entry.end if trans.entry.end > self.entry.end else self.entry.end
        # print(f"update_boundaries:  {self.entry.start}-{self.entry.end}")

    def print_gene(self):
        print(self.entry)
        for key, trans in self.transcripts.items():
            trans.print_transcript()



class Lifton_TRANS:
    def __init__(self, ref_trans_id, ref_gene_id, gene_id, copy_num, gffutil_entry_trans, ref_trans_attrs):
        ###########################
        # Assigning the reference transcripts & attributes
        ###########################
        if int(copy_num) > 0:
            trans_id = f"{ref_trans_id}_{copy_num}" 
        else:
            trans_id = ref_trans_id
        self.entry = gffutil_entry_trans
        self.entry.attributes = ref_trans_attrs
        self.entry.id = trans_id
        self.entry.attributes["ID"] = [trans_id]
        self.entry.source = "LiftOn"
        self.entry.featuretype = "transcript"
        self.exons = []
        self.exon_dic = {}

        ###########################
        # Get reference ID for the gene & transcript.
        ###########################
        self.ref_gene_id = ref_gene_id
        self.ref_tran_id = ref_trans_id

        ###########################
        # Update gene copy number dictionary
        ###########################
        # self.__update_trans_copy(self.ref_gene_id)
        self.entry.attributes['Parent'] = [gene_id]
        self.entry.attributes['transcript_id'] = [self.entry.id]


    # def __update_trans_copy(self, ref_gene_id, gene_copy_num_dict, trans_copy_num_dict):
    #     if ref_gene_id in gene_copy_num_dict.keys() and gene_copy_num_dict[ref_gene_id] > 0:
    #         trans_copy_num_dict[self.ref_tran_id] = gene_copy_num_dict[ref_gene_id]
    #         new_trans_id = self.ref_tran_id + '_' + str(trans_copy_num_dict[self.ref_tran_id])
    #         self.__update_trans_id(new_trans_id)
    #         self.__update_trans_extra_copy_attrs(trans_copy_num_dict)
    #     else:
    #         gene_copy_num_dict[ref_gene_id] = 0
    #         trans_copy_num_dict[self.ref_tran_id] = 0
    #         new_trans_id = self.ref_tran_id
    #         self.__update_trans_id(new_trans_id)

    # def __update_trans_id(self, trans_id):
    #     self.entry.id = trans_id
    #     self.entry.attributes['ID'] = [trans_id]

    # def __update_trans_extra_copy_attrs(self, trans_copy_num_dict):
    #     self.entry.attributes["extra_copy_number"] = [str(trans_copy_num_dict[self.ref_tran_id])]







    def add_lifton_status_attrs(self, lifton_status):
        self.entry.attributes["protein_identity"] = [f"{lifton_status.lifton:.3f}"]
        self.entry.attributes["status"] = [lifton_status.annotation]

    def add_exon(self, gffutil_entry_exon):
        attributes = {}
        attributes['Parent'] = [self.entry.id]
        gffutil_entry_exon.attributes = attributes
        Lifton_exon = Lifton_EXON(gffutil_entry_exon)
        lifton_utils.custom_bisect_insert(self.exons, Lifton_exon)

    def add_cds(self, gffutil_entry_cds):
        for exon in self.exons:
            if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (gffutil_entry_cds.start, gffutil_entry_cds.end)):
                attributes = {}
                attributes['Parent'] = [self.entry.id]
                gffutil_entry_cds.attributes = attributes
                exon.add_cds(gffutil_entry_cds)

    def update_gffutil_entry_trans(self, gffutil_entry_trans):
        for key, atr in gffutil_entry_trans.attributes.items():
            self.entry.attributes[key] = atr

    def update_cds_list(self, cds_list):
        idx_exon_itr = 0
        new_exons = []

        # Update frame in CDS first! => CDSs are in the protein translation order
        accum_cds_length = 0
        for cds_idx, cds in enumerate(cds_list):
            cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
            accum_cds_length = cds.entry.end - cds.entry.start + 1

        # Reverse CDS list if the strand is "-" => CDSs are in small to large order
        if self.entry.strand == "-":
            cds_list.reverse()

        ########################
        # Case 1: only 1 CDS 
        ########################
        if len(cds_list) == 1:
            only_cds = cds_list[0]
            while idx_exon_itr < len(self.exons):
                exon = self.exons[idx_exon_itr]
                idx_exon_itr += 1
                # print(f"Checking overlapping: {exon.entry.start}-{exon.entry.end}; {only_cds.entry.start}-{only_cds.entry.end}")
                if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (only_cds.entry.start, only_cds.entry.end)):

                    if exon.entry.start >= only_cds.entry.start:
                        exon.entry.start = only_cds.entry.start
                    
                    elif exon.entry.end <= only_cds.entry.end:
                        exon.entry.end = only_cds.entry.end
                    
                    exon.add_lifton_cds(only_cds)
                    new_exons.append(exon)
                else:
                    exon.add_lifton_cds(None)
                    new_exons.append(exon)

        ########################
        # Case 2: only 1 exon, and >1 CDSs
        ########################
        elif len(self.exons) == 1:
            first_exon = self.exons[0]
            for cds_idx in range(len(cds_list)):
                exon = copy.deepcopy(first_exon)
                
                cds = cds_list[cds_idx]
                if cds_idx == 0:
                    if exon.entry.start >= cds.entry.start:
                        exon.entry.start = cds.entry.start
                    exon.entry.end = cds.entry.end
                elif cds_idx == len(cds_list)-1:
                    if exon.entry.end <= cds.entry.end:
                        exon.entry.end = cds.entry.end
                    exon.entry.start = cds.entry.start
                else:
                    exon.entry.end = cds.entry.end
                    exon.entry.start = cds.entry.start
                exon.add_lifton_cds(cds)
                new_exons.append(exon)

        ########################
        # Case 3: multiple CDSs and exons
        ########################
        elif len(cds_list) > 1 and len(self.exons) > 1:
            cds_idx = 0
            exon_idx = 0
            cds = cds_list[cds_idx]
            exon = self.exons[exon_idx]
            ################################################
            # Step 1: CDSs or exons are smaller & no overlapping 
            #      => finding the first overlapping CDS and exons
            #      => processing all exons / CDS ahead of the first overlapping
            ################################################
            init_head_order = None
            if exon.entry.start >= cds.entry.end:
                init_head_order = 0
            elif exon.entry.end <= cds.entry.start:
                init_head_order = 1
                
            # print(">> Step 1: find the first overlapping exon & CDS")
            while not lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (cds.entry.start, cds.entry.end)):
                # print(f"cds_idx: {cds_idx}; {cds.entry.start}-{cds.entry.end} (len: {len(cds_list)});  exon_idx: {exon_idx}; {exon.entry.start}-{exon.entry.end} (len: {len(self.exons)})")
                if exon.entry.start >= cds.entry.end:
                    # print("Order: ", 0)
                    if init_head_order != 0:
                        # print(">> Change head break! 0")
                        break

                    # |ccccc| |eeeee|
                    # print("|ccccc| |eeeee|:")
                    # Step 1: create a new exon using CDS boundary
                    new_exon = copy.deepcopy(exon)
                    new_exon.update_exon_info(cds.entry.start, cds.entry.end)
                    # print(">> New exon!!")
                    # new_exon.print_exon()
                    
                    # Step 2: add the CDS into the exon
                    new_exon.add_lifton_cds(cds)
                    new_exons.append(new_exon)

                    # Step 3: push the first CDS to the next CDS
                    cds_idx += 1
                    cds = cds_list[cds_idx]
                elif exon.entry.end <= cds.entry.start:
                    # print("Order: ", 1)
                    if init_head_order != 1:
                        # print(">> Change head break! 1")
                        break

                    # |eeeee| |ccccc|
                    # print("|eeeee| |ccccc|:")
                    # Step 1: create a new exon using exon boundary
                    new_exon = copy.deepcopy(exon)
                    new_exon.update_exon_info(exon.entry.start, exon.entry.end)
                    # print(">> New exon!!")
                    # new_exon.print_exon()
                    new_exons.append(new_exon)

                    # Step 2: push the first exon to the next exon
                    exon_idx += 1
                    exon = self.exons[exon_idx]

            ################################################
            # Step 2: parse the first overlapping exon & CDS
            ################################################
            if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (cds.entry.start, cds.entry.end)):
                # print(">> Step 2: parse the first overlapping exon & CDS")
                # print(f"cds_idx: {cds_idx}; {cds.entry.start}-{cds.entry.end} (len: {len(cds_list)});  exon_idx: {exon_idx}; {exon.entry.start}-{exon.entry.end} (len: {len(self.exons)})")

                new_exon = copy.deepcopy(exon)
                if exon.entry.start > cds.entry.start:
                    new_exon.entry.start = cds.entry.start
                if exon.entry.end is not cds.entry.end:
                    new_exon.entry.end = cds.entry.end
                
                new_exon.add_lifton_cds(cds)
                new_exons.append(new_exon)
                # Move to the next CDS if it is not the last CDS
                if cds_idx < len(cds_list)-1:
                    cds_idx += 1

            ################################################
            # Step 3: parse the inner exons and CDSs
            ################################################
            # Handle the CDSs in the middle
            # Leave the last CDS not processed.
            # print(">> Step 3: parse the inner exons and CDSs")
            while cds_idx < len(cds_list)-1:
                cds = cds_list[cds_idx]
                # Step 1. Create a new exon
                # Step 2. Add the CDS in the new exon
                new_exon = copy.deepcopy(exon)
                new_exon.update_exon_info(cds.entry.start, cds.entry.end)
                # print(">> New exon!!")
                # new_exon.print_exon()

                new_exon.add_lifton_cds(cds)
                new_exons.append(new_exon)
                cds_idx += 1
            

            ################################################
            # Step 4: parse the last CDS & its overlapping exon
            ################################################
            cds = cds_list[cds_idx]
            exon = self.exons[exon_idx]
            
            # print(">> Step 4: parse the last CDS & its overlapping exon")
            while exon_idx < len(self.exons):
                exon = self.exons[exon_idx]
                # print(f"cds_idx: {cds_idx}; {cds.entry.start}-{cds.entry.end} (len: {len(cds_list)});  exon_idx: {exon_idx}; {exon.entry.start}-{exon.entry.end} (len: {len(self.exons)})")

                # Step 1. Create a new exon, and add cds into the exon
                new_exon = copy.deepcopy(exon)
                if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (cds.entry.start, cds.entry.end)):
                    # Step 2. Exon & CDS overlapping. If the exon is further then 
                    if new_exon.entry.end < cds.entry.end:
                        new_exon.entry.end = cds.entry.end
                    if new_exon.entry.start is not cds.entry.start:
                        new_exon.entry.start = cds.entry.start

                    # print(">>>>>>>>>>> segments_overlap")
                    # new_exon.print_exon()
                    # cds.print_cds()

                    new_exon.add_lifton_cds(cds)
                    new_exons.append(new_exon)

                else:
                    # Step 2. No overlapping. If the exon is further then 
                    # print(">>>>>>>>>>> segments do not overlap")
                    if new_exon.entry.start >= cds.entry.end:
                        # |ccccc| |eeeee|
                        # All the subsequent exons are after the last CDS => create new exon & CDS for all of them!

                        # print("#############")
                        # new_exon.print_exon()
                        new_exon.add_lifton_cds(None)
                        new_exons.append(new_exon)
                    elif exon_idx == len(self.exons)-1 and cds.entry.start >= new_exon.entry.end:
                        # |eeeee| |ccccc| 
                        # All the last exons are after before the last CDS => create the last CDS!
                        new_exon.update_exon_info(cds.entry.start, cds.entry.end)
                        # print(">> !!New exon!!")
                        new_exon.add_lifton_cds(cds)
                        new_exons.append(new_exon)
                # new_exon.print_exon()
                exon_idx += 1

        self.exons = new_exons
        self.update_boundaries()


    def fix_truncated_protein(self, fai, ref_protein_seq, ref_trans_seq, lifton_status):
        # Need to output which type of mutation it is.
        ################################
        # Step 1: Iterate through the children and chain the DNA sequence
        ################################
        coding_seq = ""
        cdss_lens = []
        trans_seq = ""
        exon_lens = []
        for exon in self.exons:
            # Chaining the exon features
            p_trans_seq = exon.entry.sequence(fai)
            p_trans_seq = Seq(p_trans_seq).upper()
            if exon.entry.strand == '-':
                trans_seq = p_trans_seq + trans_seq
                exon_lens.insert(0, exon.entry.end - exon.entry.start + 1)
            elif exon.entry.strand == '+':
                trans_seq = trans_seq + p_trans_seq
                exon_lens.append(exon.entry.end - exon.entry.start + 1)

            if exon.cds is not None:
                # Chaining the CDS features
                p_seq = exon.cds.entry.sequence(fai)
                p_seq = Seq(p_seq).upper()
                if exon.cds.entry.strand == '-':
                    coding_seq = p_seq + coding_seq
                    cdss_lens.insert(0, exon.cds.entry.end - exon.cds.entry.start + 1)
                elif exon.cds.entry.strand == '+':
                    coding_seq = coding_seq + p_seq
                    cdss_lens.append(exon.cds.entry.end - exon.cds.entry.start + 1)

        ################################
        # Step 2: Translate the DNA sequence & get the reference protein sequence.
        ################################
        protein_seq = coding_seq.translate()
        peps = protein_seq.split("*")

        lifton_aa_aln = align.protein_align(protein_seq, ref_protein_seq)
        # Update lifton sequence identity
        lifton_status.lifton = max(lifton_status.lifton, lifton_aa_aln.identity)

        lifton_tran_aln = align.trans_align(coding_seq, ref_trans_seq)
        variants.find_variants(lifton_tran_aln, lifton_aa_aln, lifton_status, peps)

        # self.entry.attributes["mutation"] 

        ORF_search = False
        for mutation in lifton_status.status:
            # identical
            # synonymous # inframe_insertion # inframe_deletion # nonsynonymous # frameshift # start_lost # stop_missing # stop_codon_gain            
            # Adding mutations in to entry.attributes
            if mutation != "identical":
                if "mutation" not in self.entry.attributes:
                    self.entry.attributes["mutation"] = [mutation]
                else:
                    self.entry.attributes["mutation"].append(mutation)
            # ORF searching
            if mutation == "stop_missing" or mutation == "stop_codon_gain" or mutation == "frameshift" or mutation == "start_lost":
                ORF_search = True

        if ORF_search:
            self.__find_orfs(trans_seq, exon_lens, ref_protein_seq, lifton_aa_aln, lifton_status)
        return lifton_tran_aln, lifton_aa_aln


    def __find_orfs(self, trans_seq, exon_lens, ref_protein_seq, lifton_aln, lifton_status):
        trans_seq = trans_seq.upper()
        # Find ORFs manually
        start_codon = "ATG"
        stop_codons = ["TAA", "TAG", "TGA"]

        orf_list = []
        max_orf_len = [0, 0, 0]
        for frame in range(3):
            orf_idx_s = 0
            for i in range(frame, len(trans_seq), 3):
                codon = str(trans_seq[i:i+3])
                if codon == start_codon:
                    orf_idx_s = i
                    orf_idx_e = i
                    orf_seq = ""
                    for j in range(i, len(trans_seq), 3):
                        codon = str(trans_seq[j:j+3])
                        orf_seq += codon
                        if codon in stop_codons:
                            orf_idx_e = j+3
                            break

                    if orf_seq and (orf_idx_s < orf_idx_e):
                        curr_orf_len = orf_idx_e-orf_idx_s+1
                        if curr_orf_len >  max_orf_len[frame]:
                            max_orf_len[frame] = curr_orf_len
                            orf = lifton_class.Lifton_ORF(orf_idx_s, orf_idx_e)
                            orf_list.append(orf)

        # Print the ORFs
        final_orf = None
        max_identity = 0
        for i, orf in enumerate(orf_list):
            orf_DNA_seq = trans_seq[orf.start:orf.end]
            orf_protein_seq = orf_DNA_seq.translate()
            
            # print(f"\tORF {i+1}: {orf_DNA_seq}")
            # print(f"\torf_protein_seq: {orf_protein_seq}")
            extracted_parasail_res, extracted_seq, reference_seq = align.parasail_align_protein_base(orf_protein_seq, str(ref_protein_seq))        

            alignment_score = extracted_parasail_res.score
            alignment_query = extracted_parasail_res.traceback.query
            alignment_comp = extracted_parasail_res.traceback.comp
            alignment_ref = extracted_parasail_res.traceback.ref

            extracted_matches, extracted_length = get_id_fraction.get_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query, 0, len(extracted_parasail_res.traceback.ref))

            extracted_identity = extracted_matches/extracted_length


            # print(f"\textracted_identity: {extracted_identity}")
            if extracted_identity > max_identity:
                max_identity = extracted_identity
                final_orf = orf
                lifton_aln = lifton_class.Lifton_Alignment(extracted_identity, None, alignment_query, alignment_comp, alignment_ref, None, None, extracted_seq, reference_seq, None)
                lifton_status.lifton = max(lifton_status.lifton, lifton_aln.identity)

        if final_orf is not None:
            self.__update_cds_boundary(final_orf)
        

    def __update_cds_boundary(self, final_orf):
        if self.entry.strand == "+":
            self.__iterate_exons_update_cds(final_orf, self.exons, "+")
        elif self.entry.strand == "-":
            self.__iterate_exons_update_cds(final_orf, self.exons[::-1], "-")            


    def __iterate_exons_update_cds(self, final_orf, exons, strand):
        accum_exon_length = 0
        accum_cds_length = 0
        for exon_idx, exon in enumerate(exons):
            curr_exon_len = exon.entry.end - exon.entry.start + 1
            # print(f"\t>> {exon.entry.seqid} {exon.entry.strand}; exon_idx: {exon_idx}: {exon.entry.start}-{exon.entry.end} (len: {len(exons)});  accum_exon_length: {accum_exon_length}; curr_exon_len: {curr_exon_len}; final_orf.start: {final_orf.start}; final_orf.end: {final_orf.end}")

            if accum_exon_length <= final_orf.start:
                if final_orf.start < accum_exon_length+curr_exon_len:
                    # Create first partial CDS
                    if exon.cds is not None:
                        if strand == "+":
                            exon.cds.entry.start = exon.entry.start + (final_orf.start - accum_exon_length)
                        elif strand == "-":
                            exon.cds.entry.end = exon.entry.end - (final_orf.start - accum_exon_length)
                    else:
                        if strand == "+":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.start + (final_orf.start - accum_exon_length), exon.entry.end)
                        elif strand == "-":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.start, exon.entry.end - (final_orf.start - accum_exon_length))
                    exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                    accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)
                    # print(f"\t\t >> {exon.entry.seqid} {exon.entry.strand}; exon_idx: {exon_idx}: {exon.entry.start}-{exon.entry.end} (len: {len(exons)});  accum_exon_length: {accum_exon_length}; curr_exon_len: {curr_exon_len}; final_orf.start: {final_orf.start}; final_orf.end: {final_orf.end}")
                    # print(f"\t\t >> exon.cds.entry.start: {exon.cds.entry.start}; exon.cds.entry.end: {exon.cds.entry.end}")
                    # print(f"\t\t>> exon.cds.entry.frame: {exon.cds.entry.frame}")
                else:
                    # No CDS should be created
                    exon.cds = None

            elif final_orf.start < accum_exon_length and accum_exon_length < final_orf.end:

                if final_orf.end <= accum_exon_length+curr_exon_len:
                    # Create the last partial CDS
                    if exon.cds is not None:
                        if strand == "+":
                            exon.cds.entry.end = exon.entry.start + (final_orf.end - accum_exon_length)-1
                        elif strand == "-":
                            exon.cds.entry.start = exon.entry.end - (final_orf.end - accum_exon_length)+1
                    else:
                        if strand == "+":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.start, exon.entry.start + (final_orf.end - accum_exon_length)-1)
                        elif strand == "-":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.end - (final_orf.end - accum_exon_length)+1, exon.entry.end)

                    exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                    accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)
                    # print(f"\t\t >> {exon.entry.seqid} {exon.entry.strand}; exon_idx: {exon_idx}: {exon.entry.start}-{exon.entry.end} (len: {len(exons)});  accum_exon_length: {accum_exon_length}; curr_exon_len: {curr_exon_len}; final_orf.start: {final_orf.start}; final_orf.end: {final_orf.end}")
                    # print(f"\t\t >> exon.cds.entry.start: {exon.cds.entry.start}; exon.cds.entry.end: {exon.cds.entry.end}")
                    # print(f"\t\t >> exon.cds.entry.frame: {exon.cds.entry.frame}")

                else:
                    # Keep the original full CDS
                    if exon.cds is None:
                        exon.add_novel_lifton_cds(exon.entry, exon.entry.start, exon.entry.end)
                    exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                    accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)

            elif final_orf.end <= accum_exon_length:
                # No CDS should be created
                exon.cds = None
            accum_exon_length += curr_exon_len

    def __get_cds_frame(self, accum_cds_length):
        return (3 - accum_cds_length%3)%3


    def write_entry(self, fw):
        # print("Inside 'write_entry'!")
        # print(self.entry)
        # print(f"(write_entry) update_boundaries, exon length: {len(self.exons)};  {self.entry.start} - {self.entry.end}")

        fw.write(str(self.entry) + "\n")
        # Write out the exons first
        for exon in self.exons:
            exon.write_entry(fw)
        # Write out the CDSs second
        for exon in self.exons:
            if exon.cds is not None:
                exon.cds.write_entry(fw)

    def update_boundaries(self):
        self.entry.start = self.exons[0].entry.start
        self.entry.end = self.exons[-1].entry.end
        # print(f"update_boundaries, exon length: {len(self.exons)};  {self.entry.start} - {self.entry.end}")

    def print_transcript(self):
        print(f"\t{self.entry}")
        for exon in self.exons:
            exon.print_exon()



class Lifton_EXON:
    def __init__(self, gffutil_entry_exon):
        gffutil_entry_exon.source = "LiftOn"
        gffutil_entry_exon.featuretype = "exon"
        self.entry = gffutil_entry_exon
        if 'extra_copy_number' in self.entry.attributes: self.entry.attributes.pop('extra_copy_number')
        self.cds = None

        # print("self.entry: ", self.entry.featuretype)

    def update_exon_info(self, start, end):
        self.cds = None
        self.entry.source = "LiftOn"
        self.entry.start = start
        self.entry.end = end
        # print(f"start: {self.entry.start} - end: {self.entry.end}")

    def add_cds(self, gffutil_entry_cds):
        Lifton_cds = Lifton_CDS(gffutil_entry_cds)
        self.cds = Lifton_cds

    def add_novel_lifton_cds(self, gffutil_entry_exon, start, end):
        gffutil_entry_cds = copy.deepcopy(gffutil_entry_exon)
        gffutil_entry_cds.featuretype = "CDS"
        gffutil_entry_cds.start = start
        gffutil_entry_cds.end = end
        Lifton_cds = Lifton_CDS(gffutil_entry_cds)
        attributes = {}
        attributes['Parent'] = self.entry.attributes['Parent']
        Lifton_cds.entry.attributes = attributes
        self.cds = Lifton_cds

    def add_lifton_cds(self, Lifton_cds):
        if Lifton_cds is not None:
            attributes = {}
            attributes['Parent'] = self.entry.attributes['Parent']
            Lifton_cds.entry.attributes = attributes
        self.cds = Lifton_cds

    def write_entry(self, fw):
        fw.write(str(self.entry) + "\n")

    def print_exon(self):
        print(f"\t\t{self.entry}")
        if self.cds != None:
            self.cds.print_cds()



class Lifton_CDS:
    def __init__(self, gffutil_entry_cds):
        gffutil_entry_cds.source = "LiftOn"
        gffutil_entry_cds.featuretype = "CDS"
        self.entry = gffutil_entry_cds
        if 'extra_copy_number' in self.entry.attributes: self.entry.attributes.pop('extra_copy_number')

    def write_entry(self, fw):
        fw.write(str(self.entry) + "\n")

    def print_cds(self):
        print(f"\t\t{self.entry}")
