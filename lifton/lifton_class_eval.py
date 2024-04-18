from lifton import lifton_class, lifton_utils, align, get_id_fraction, variants
import copy, os
from Bio.Seq import Seq


class Lifton_ORF_EVAL:
    def __init__(self, start, end):
        self.start = start
        self.end = end


class Lifton_GENE_EVAL:
    def __init__(self, ref_gene_id, gffutil_entry_gene):
        ###########################
        # Assigning the reference gene & attributes
        ###########################
        self.entry = gffutil_entry_gene
        self.transcripts = {}
        self.ref_gene_id = ref_gene_id

    def add_transcript(self, ref_trans_id, gffutil_entry_trans):
        Lifton_trans = Lifton_TRANS_EVAL(ref_trans_id, self.ref_gene_id, self.entry.id, gffutil_entry_trans)
        self.transcripts[Lifton_trans.entry.id] = Lifton_trans
        return Lifton_trans


    def add_exon(self, trans_id, gffutil_entry_exon):
        self.transcripts[trans_id].add_exon(gffutil_entry_exon)

    def add_cds(self, trans_id, gffutil_entry_cds):
        self.transcripts[trans_id].add_cds(gffutil_entry_cds)

    def add_feature(self, gffutil_entry_trans):
        Lifton_feature = LiftOn_FEATURE_EVAL(self.entry.id, gffutil_entry_trans)
        self.transcripts[Lifton_feature.entry.id] = Lifton_feature
        return Lifton_feature
    
    def fix_truncated_protein(self, trans_id, ref_trans_id, fai, ref_proteins, fai_trans, lifton_status):
        ref_protein_seq = ref_proteins[ref_trans_id]
        ref_trans_seq = fai_trans[ref_trans_id]
        if trans_id not in self.transcripts.keys():
            return None, False
        lifton_aln, good_trans = self.transcripts[trans_id].fix_truncated_protein(fai, ref_protein_seq, ref_trans_seq, lifton_status)
        return lifton_aln, good_trans
               

class LiftOn_FEATURE_EVAL:
    def __init__(self, parent_id, gffutil_entry_feature):
        self.entry = gffutil_entry_feature
        self.features = {}
        self.entry.attributes["Parent"] = [parent_id]

    def add_feature(self, gffutil_entry_trans):
        Lifton_feature = LiftOn_FEATURE_EVAL(self.entry.id, gffutil_entry_trans)
        self.features[Lifton_feature.entry.id] = Lifton_feature
        return Lifton_feature
    

class Lifton_TRANS_EVAL:
    def __init__(self, ref_trans_id, ref_gene_id, gene_id, gffutil_entry_trans):

        ###########################
        # Assigning the reference transcripts & attributes
        ###########################
        self.entry = gffutil_entry_trans
        self.exons = []
        self.exon_dic = {}

        ###########################
        # Get reference ID for the gene & transcript.
        ###########################
        self.ref_gene_id = ref_gene_id
        self.ref_tran_id = ref_trans_id

    def add_exon(self, gffutil_entry_exon):
        attributes = {}
        attributes['Parent'] = [self.entry.id]
        gffutil_entry_exon.attributes = attributes
        Lifton_exon = Lifton_EXON_EVAL(gffutil_entry_exon)
        lifton_utils.custom_bisect_insert(self.exons, Lifton_exon)

    def add_cds(self, gffutil_entry_cds):
        for exon in self.exons:
            _, ovp = lifton_utils.segments_overlap_length((exon.entry.start, exon.entry.end), (gffutil_entry_cds.start, gffutil_entry_cds.end))
            if ovp:
                attributes = {}
                attributes['Parent'] = [self.entry.id]
                gffutil_entry_cds.attributes = attributes
                exon.add_cds(gffutil_entry_cds)

    def update_gffutil_entry_trans(self, gffutil_entry_trans):
        for key, atr in gffutil_entry_trans.attributes.items():
            self.entry.attributes[key] = atr


    def get_coding_seq(self, fai):
        coding_seq = ""
        cdss_lens = []
        cds_children = []
        for exon in self.exons:
            if exon.cds is not None:
                cds_children.append(copy.deepcopy(exon.cds.entry))
                # Chaining the CDS features
                p_seq = exon.cds.entry.sequence(fai)
                if exon.cds.entry.strand == '-':
                    coding_seq = p_seq + coding_seq
                    cdss_lens.insert(0, exon.cds.entry.end - exon.cds.entry.start + 1)
                elif exon.cds.entry.strand == '+':
                    coding_seq = coding_seq + p_seq
                    cdss_lens.append(exon.cds.entry.end - exon.cds.entry.start + 1)
        return coding_seq, cds_children, cdss_lens

    def get_coding_trans_seq(self, fai):
        trans_seq = ""
        coding_seq = ""
        accum_cds_length = 0
        for exon in self.exons:
            # Chaining the exon features
            p_trans_seq = exon.entry.sequence(fai)
            p_trans_seq = Seq(p_trans_seq).upper()
            if exon.entry.strand == '-':
                trans_seq = p_trans_seq + trans_seq
            elif exon.entry.strand == '+':
                trans_seq = trans_seq + p_trans_seq
            if exon.cds is not None:
                # Updating CDS frame
                exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                accum_cds_length = exon.cds.entry.end - exon.cds.entry.start + 1
                # Chaining the CDS features
                p_seq = exon.cds.entry.sequence(fai)
                if exon.cds.entry.strand == '-':
                    coding_seq = p_seq + coding_seq
                elif exon.cds.entry.strand == '+':
                    coding_seq = coding_seq + p_seq
        if trans_seq != None:
            trans_seq = str(trans_seq).upper()
        if coding_seq != None:
            coding_seq = str(coding_seq).upper()
        return coding_seq, trans_seq

    def translate_coding_seq(self, coding_seq):
        protein_seq = None
        if coding_seq != "":
            protein_seq = str(Seq(coding_seq).translate())
        return protein_seq

    def align_coding_seq(self, protein_seq, ref_protein_seq, lifton_status):
        if protein_seq == "" or protein_seq == None:
            lifton_aa_aln = None
            peps = None
        else:
            peps = protein_seq.split("*")
            lifton_aa_aln = align.protein_align(str(protein_seq), str(ref_protein_seq))
            # Update lifton sequence identity
            lifton_status.lifton_aa = max(lifton_status.lifton_aa, lifton_aa_aln.identity)
        return lifton_aa_aln, peps

    def align_trans_seq(self, trans_seq, ref_trans_seq, lifton_status):
        if trans_seq == "" or trans_seq == None:
            lifton_tran_aln = None
        else:
            lifton_tran_aln = align.trans_align(trans_seq, ref_trans_seq)
        lifton_status.lifton_dna = lifton_tran_aln.identity
        return lifton_tran_aln

    def fix_truncated_protein(self, fai, ref_protein_seq, ref_trans_seq, lifton_status):
        # Getting translated sequences (frame will be updated)
        coding_seq, trans_seq = self.get_coding_trans_seq(fai)
        protein_seq = self.translate_coding_seq(coding_seq)
        # Aligning the LiftOn protein & DNA sequences
        lifton_aa_aln, peps = self.align_coding_seq(protein_seq, ref_protein_seq, lifton_status)
        lifton_tran_aln = self.align_trans_seq(str(trans_seq), str(ref_trans_seq), lifton_status)
        variants.find_variants(lifton_tran_aln, lifton_aa_aln, lifton_status, peps)
        ORF_search = False
        for mutation in lifton_status.status:
            # identical # synonymous 
            # (1) inframe_insertion # (2) inframe_deletion # (3) nonsynonymous 
            # (4) frameshift # (5) start_lost # (6) stop_missing # (7) stop_codon_gain
            # Adding mutations in to entry.attributes
            if mutation != "identical":
                if "mutation" not in self.entry.attributes:
                    self.entry.attributes["mutation"] = [mutation]
                else:
                    self.entry.attributes["mutation"].append(mutation)
            # ORF searching for these four types of mutations
            # frameshift_orf_threshold = 0.8
            # and lifton_aa_aln.identity < frameshift_orf_threshold)
            if mutation == "stop_missing" or mutation == "stop_codon_gain" or mutation == "frameshift"  or mutation == "start_lost":
                ORF_search = True
        # if ORF_search:
        #     self.__find_orfs(trans_seq, ref_protein_seq, lifton_aa_aln, lifton_status)
        return lifton_tran_aln, lifton_aa_aln

    def align_trans_dna(self, fai, ref_trans_seq, lifton_status):
        _, trans_seq = self.get_coding_trans_seq(fai)
        lifton_tran_aln = self.align_trans_seq(trans_seq, ref_trans_seq, lifton_status)
        return lifton_tran_aln

    def __find_orfs(self, trans_seq, ref_protein_seq, lifton_aln, lifton_status):
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
            print(f"\tORF {i+1}: {orf.start}-{orf.end}")

            orf_DNA_seq = trans_seq[orf.start:orf.end]
            orf_protein_seq = orf_DNA_seq.translate()
            
            # print(f"\tORF {i+1}: {orf_DNA_seq}")
            # print(f"\torf_protein_seq: {orf_protein_seq}")
            extracted_parasail_res = align.parasail_align_protein_base(orf_protein_seq,ref_protein_seq)

            alignment_score = extracted_parasail_res.score
            alignment_query = extracted_parasail_res.traceback.query
            alignment_comp = extracted_parasail_res.traceback.comp
            alignment_ref = extracted_parasail_res.traceback.ref

            extracted_matches, extracted_length = get_id_fraction.get_AA_id_fraction(extracted_parasail_res.traceback.ref, extracted_parasail_res.traceback.query)

            extracted_identity = extracted_matches/extracted_length

            # print(f"\textracted_identity: {extracted_identity}")
            if extracted_identity > max_identity:
                max_identity = extracted_identity
                final_orf = orf
                lifton_aln = lifton_class.Lifton_Alignment(extracted_identity, None, alignment_query, alignment_comp, alignment_ref, None, None, orf_protein_seq, ref_protein_seq, None)
                lifton_status.lifton_aa = lifton_aln.identity
                # max(lifton_status.lifton_aa, lifton_aln.identity)

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


class Lifton_EXON_EVAL:
    def __init__(self, gffutil_entry_exon):
        self.entry = gffutil_entry_exon
        self.cds = None

    def add_cds(self, gffutil_entry_cds):
        Lifton_cds = Lifton_CDS_EVAL(gffutil_entry_cds)
        self.cds = Lifton_cds

    def add_lifton_cds(self, Lifton_cds):
        if Lifton_cds is not None:
            attributes = {}
            attributes['Parent'] = self.entry.attributes['Parent']
            Lifton_cds.entry.attributes = attributes
        self.cds = Lifton_cds

    def add_novel_lifton_cds(self, gffutil_entry_exon, start, end):
        gffutil_entry_cds = copy.deepcopy(gffutil_entry_exon)
        gffutil_entry_cds.featuretype = "CDS"
        gffutil_entry_cds.start = start
        gffutil_entry_cds.end = end
        Lifton_cds = Lifton_CDS_EVAL(gffutil_entry_cds)
        attributes = {}
        attributes['Parent'] = self.entry.attributes['Parent']
        Lifton_cds.entry.attributes = attributes
        self.cds = Lifton_cds


class Lifton_CDS_EVAL:
    def __init__(self, gffutil_entry_cds):
        self.entry = gffutil_entry_cds