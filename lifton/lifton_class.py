from lifton import align, lifton_utils, lifton_class, get_id_fraction, variants, logger
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
        self.lifton_dna = 0
        self.lifton_aa = 0
        self.eval_dna = 0
        self.eval_aa = 0
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
        self.is_protein_coding = False
        self.is_non_coding = False
        self.children = set()


class Lifton_GENE:
    def __init__(self, ref_gene_id, gffutil_entry_gene, ref_gene_attrs, tree_dict, ref_features_dict, args, tmp = False):
        ###########################
        # Assigning the reference gene & attributes
        ###########################
        self.entry = gffutil_entry_gene
        self.entry.source = "LiftOn"
        self.is_protein_coding = False
        self.is_non_coding = False
        self.transcripts = {}
        self.ref_gene_id = ref_gene_id
        self.copy_num = self.__get_gene_copy(ref_features_dict)
        self.tmp = tmp
        self.entry.attributes = ref_gene_attrs
        self.entry.attributes["ID"] = self.ref_gene_id + "_" + str(self.copy_num) if self.copy_num > 0 else self.ref_gene_id
        if self.copy_num > 0:
            self.entry.attributes["extra_copy_number"] = [str(self.copy_num)]        
        self.__update_gene_copy(ref_features_dict)
        self.entry.id = self.entry.attributes["ID"][0]
        gene_interval = Interval(self.entry.start, self.entry.end, self.entry.id)
        if self.entry.seqid not in tree_dict.keys():
            tree_dict[self.entry.seqid] = IntervalTree()
        tree_dict[self.entry.seqid].add(gene_interval)
        # Decide if its type
        gene_type_key = ""
        if args.annotation_database.upper() == "REFSEQ":
            gene_type_key = "gene_biotype"
        elif args.annotation_database.upper() == "GENCODE" or args.annotation_database.upper() == "ENSEMBL" or args.annotation_database.upper() == "CHESS":
            gene_type_key = "gene_type"
        if gene_type_key in self.entry.attributes.keys():
            if self.entry.attributes[gene_type_key][0] == "protein_coding":
                self.is_protein_coding = True
            elif (self.entry.attributes[gene_type_key][0] == "lncRNA" or self.entry.attributes[gene_type_key][0] == "ncRNA"):
                self.is_non_coding = True
        
    def __get_gene_copy(self, ref_features_dict):
        if 'extra_copy_number' in self.entry.attributes:
            return int(self.entry.attributes['extra_copy_number'][0])
        else:
            if self.ref_gene_id in ref_features_dict.keys():
                return ref_features_dict[self.ref_gene_id].copy_num
            else:
                return 0

    def __update_gene_copy(self, ref_features_dict):
        if self.ref_gene_id in ref_features_dict.keys():
            ref_features_dict[self.ref_gene_id].copy_num += 1

    def update_gene_info(self, chromosome, start, end):
        self.entry.seqid = chromosome
        self.entry.start = start
        self.entry.end = end

    def add_miniprot_transcript(self, ref_trans_id, gffutil_entry_trans, ref_trans_attrs, ref_features_dict):
        Lifton_trans = Lifton_TRANS(ref_trans_id, self.ref_gene_id, self.entry.id, self.copy_num, gffutil_entry_trans, ref_trans_attrs)
        self.transcripts[Lifton_trans.entry.id] = Lifton_trans
        return Lifton_trans
    
    def update_trans_info(self, trans_id, chromosome, start, end):
        self.transcripts[trans_id].entry.seqid = chromosome
        self.transcripts[trans_id].entry.start = start
        self.transcripts[trans_id].entry.end = end

    def add_transcript(self, ref_trans_id, gffutil_entry_trans, ref_trans_attrs):
        Lifton_trans = Lifton_TRANS(ref_trans_id, self.ref_gene_id, self.entry.id, self.copy_num, gffutil_entry_trans, ref_trans_attrs)
        self.transcripts[Lifton_trans.entry.id] = Lifton_trans
        return Lifton_trans

    def add_feature(self, gffutil_entry_trans):
        Lifton_feature = LiftOn_FEATURE(self.entry.id, gffutil_entry_trans, self.copy_num)
        self.transcripts[Lifton_feature.entry.id] = Lifton_feature
        return Lifton_feature

    def add_exon(self, trans_id, gffutil_entry_exon):
        self.transcripts[trans_id].add_exon(gffutil_entry_exon)

    def add_cds(self, trans_id, gffutil_entry_cds):
        self.transcripts[trans_id].add_cds(gffutil_entry_cds)
                            
    def orf_search_protein(self, trans_id, ref_trans_id, fai, ref_proteins, fai_trans, lifton_status, eval_only=False, eval_liftoff_chm13=False):
        if trans_id not in self.transcripts.keys():
            return None, False
        if not eval_liftoff_chm13:
            ref_protein_seq = str(ref_proteins[ref_trans_id]) if ref_trans_id in ref_proteins.keys() else None
            ref_trans_seq = str(fai_trans[ref_trans_id]) if ref_trans_id in fai_trans.keys() else None
        else:
            ref_protein_seq = str(ref_proteins["rna-"+ref_trans_id]) if "rna-"+ref_trans_id in ref_proteins.keys() else None
            ref_trans_seq = str(fai_trans["rna-"+ref_trans_id]) if "rna-"+ref_trans_id in fai_trans.keys() else None
        lifton_aln, good_trans = self.transcripts[trans_id].orf_search_protein(fai, ref_protein_seq, ref_trans_seq, lifton_status, is_non_coding=self.is_non_coding, eval_only=eval_only)
        return lifton_aln, good_trans

    def update_cds_list(self, trans_id, cds_list):
        self.transcripts[trans_id].update_cds_list(cds_list)
        self.update_boundaries()

    def add_lifton_gene_status_attrs(self, source):
        self.entry.attributes["source"] = [source]

    def add_lifton_trans_status_attrs(self, trans_id, lifton_status):
        self.transcripts[trans_id].add_lifton_trans_status_attrs(lifton_status)

    def write_entry(self, fw, transcripts_stats_dict):
        if not self.tmp:
            fw.write(str(self.entry) + "\n")
        for key, trans in self.transcripts.items():
            trans.write_entry(fw)
            TYPE = ""
            if self.is_protein_coding and trans.entry.featuretype == "mRNA":
                TYPE = "coding"
            elif self.is_non_coding and (trans.entry.featuretype == "ncRNA" or trans.entry.featuretype == "nc_RNA" or trans.entry.featuretype == "lncRNA" or trans.entry.featuretype == "lnc_RNA"):
                TYPE = "non-coding"
            else:
                TYPE = "other"
            if trans.ref_tran_id is None:
                continue
            if not trans.ref_tran_id in transcripts_stats_dict[TYPE].keys():
                transcripts_stats_dict[TYPE][trans.ref_tran_id] = 1
            else:
                transcripts_stats_dict[TYPE][trans.ref_tran_id] += 1

    def update_boundaries(self):        
        for key, trans in self.transcripts.items():
            self.entry.start = trans.entry.start if trans.entry.start < self.entry.start else self.entry.start
            self.entry.end = trans.entry.end if trans.entry.end > self.entry.end else self.entry.end

    def print_gene(self):
        print(self.entry)
        for key, trans in self.transcripts.items():
            trans.print_transcript()


class LiftOn_FEATURE:
    def __init__(self, parent_id, gffutil_entry_feature, copy_num):
        self.entry = gffutil_entry_feature
        self.entry.source = "LiftOn"
        self.copy_num = copy_num
        self.features = {}
        self.entry.attributes["Parent"] = [parent_id]
        self.ref_tran_id = None
        # self.ref_tran_id = self.entry.id
        if int(copy_num) > 0:
            feature_id_base = lifton_utils.get_ID_base(self.entry.id)
            self.entry.id = f"{feature_id_base}_{copy_num}"
            self.entry.attributes["ID"] = [self.entry.id]

    def add_feature(self, gffutil_entry_trans):
        Lifton_feature = LiftOn_FEATURE(self.entry.id, gffutil_entry_trans, self.copy_num)
        self.features[Lifton_feature.entry.id] = Lifton_feature
        return Lifton_feature

    def write_entry(self, fw):
        fw.write(str(self.entry) + "\n")
        for key, feature in self.features.items():
            feature.write_entry(fw)

    def print_feature(self):
        print(self.entry)
        for key, feature in self.features.items():
            feature.print_feature()


class Lifton_TRANS:
    def __init__(self, ref_trans_id, ref_gene_id, gene_id, copy_num, gffutil_entry_trans, ref_trans_attrs):
        # Assigning the reference transcripts & attributes
        if int(copy_num) > 0:
            trans_id = f"{ref_trans_id}_{copy_num}" 
        else:
            trans_id = ref_trans_id
        self.entry = gffutil_entry_trans
        self.entry.attributes = ref_trans_attrs
        self.entry.id = trans_id
        self.entry.attributes["ID"] = [trans_id]
        self.entry.source = "LiftOn"
        self.exons = []
        self.exon_dic = {}
        # Get reference ID for the gene & transcript.
        self.ref_gene_id = ref_gene_id
        self.ref_tran_id = ref_trans_id
        if 'Parent' in self.entry.attributes:
            self.entry.attributes['Parent'] = [gene_id]
        if 'transcript_id' in self.entry.attributes:
            self.entry.attributes['transcript_id'] = [self.entry.id]

    def add_lifton_trans_status_attrs(self, lifton_status):
        self.entry.attributes["protein_identity"] = [f"{lifton_status.lifton_aa:.3f}"]
        self.entry.attributes["dna_identity"] = [f"{lifton_status.lifton_dna:.3f}"]
        self.entry.attributes["status"] = [lifton_status.annotation]

    def add_exon(self, gffutil_entry_exon):
        gffutil_entry_exon.attributes['Parent'] = [self.entry.id]
        Lifton_exon = Lifton_EXON(gffutil_entry_exon)
        lifton_utils.custom_bisect_insert(self.exons, Lifton_exon)

    def add_cds(self, gffutil_entry_cds):
        for exon in self.exons:
            _, ovp = lifton_utils.segments_overlap_length((exon.entry.start, exon.entry.end), (gffutil_entry_cds.start, gffutil_entry_cds.end))
            if ovp:
                gffutil_entry_cds.attributes['Parent'] = [self.entry.id]
                exon.add_cds(gffutil_entry_cds)

    def update_gffutil_entry_trans(self, gffutil_entry_trans):
        for key, atr in gffutil_entry_trans.attributes.items():
            self.entry.attributes[key] = atr

    def mv_exon_idx(self, exon_idx):
        if exon_idx < len(self.exons)-1:
            exon_idx += 1
        return exon_idx

    def update_cds_list(self, cds_list):
        idx_exon_itr = 0
        new_exons = []
        # Reverse CDS list if the strand is "-" => CDSs are in small to large order
        if self.entry.strand == "-":
            cds_list.reverse()
        ########################
        # Case 1: only 1 CDS 
        ########################
        if len(cds_list) == 1:
            only_cds = cds_list[0]
            processed_ovp_exons = False
            ovp_exons = []
            while idx_exon_itr < len(self.exons):
                exon = self.exons[idx_exon_itr]
                _, ovp = lifton_utils.segments_overlap_length((exon.entry.start, exon.entry.end), (only_cds.entry.start, only_cds.entry.end))

                # |eeeeee| |ccccc|
                if exon.entry.end < only_cds.entry.start:
                    new_exons.append(exon)
                elif ovp:
                    cp_exon = copy.deepcopy(exon)
                    ovp_exons.append(cp_exon)
                # |cccc|  |eeee|
                elif exon.entry.start > only_cds.entry.end: 
                    processed_ovp_exons = True
                    merged_exon = copy.deepcopy(exon)
                    if len(ovp_exons) == 0:
                        cp_exon = copy.deepcopy(exon)
                        ovp_exons.append(cp_exon)
                    merged_exon.entry.start = ovp_exons[0].entry.start
                    merged_exon.entry.end = ovp_exons[-1].entry.end
                    merged_exon.add_lifton_cds(only_cds)
                    new_exons.append(merged_exon)
                    new_exons.append(exon)
                idx_exon_itr += 1
            if not processed_ovp_exons:
                merged_exon = copy.deepcopy(exon)
                merged_exon.entry.start = ovp_exons[0].entry.start
                merged_exon.entry.end = ovp_exons[-1].entry.end
                merged_exon.add_lifton_cds(only_cds)
                new_exons.append(merged_exon)
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
            init_head_order = False
            # |ccccc||eeeeee|
            if exon.entry.start >= cds.entry.end:
                init_head_order = False
            # |eeeeee| |ccccc|
            elif exon.entry.start < cds.entry.start:
                init_head_order = True
            if init_head_order:
                # Find the first overlapping exon with CDS
                while exon_idx < len(self.exons) and cds_idx < len(cds_list):
                    exon_start, exon_end = self.exons[exon_idx].entry.start, self.exons[exon_idx].entry.end
                    cds_start, cds_end = cds_list[cds_idx].entry.start, cds_list[cds_idx].entry.end
                    _, ovp = lifton_utils.segments_overlap_length((exon_start, exon_end), (cds_start, cds_end))
                    # |eeeeee| |ccccc|
                    if exon_end < cds_start:
                        # Create a new exon using exon boundary
                        new_exon = copy.deepcopy(self.exons[exon_idx])
                        new_exon.update_exon_info(exon_start, exon_end)
                        new_exons.append(new_exon)
                        exon_idx += 1
                    # |eeee|--|cccc|
                    elif ovp:
                        # Stop if exon and CDS overlap
                        # Create a new exon using exon boundary
                        new_exon = copy.deepcopy(self.exons[exon_idx])
                        new_exon.update_exon_info(min(exon_start, cds_start), cds_end)
                        new_exon.add_lifton_cds(cds_list[0])
                        new_exons.append(new_exon)
                        exon_idx += 1
                        cds_idx += 1
                        break
                    # |cccc|  |eeee|
                    elif exon_start > cds_end: 
                        # Create a new exon using CDS boundary
                        new_exon = copy.deepcopy(self.exons[exon_idx])
                        new_exon.update_exon_info(cds_start, cds_end)
                        new_exon.add_lifton_cds(cds_list[0])
                        new_exons.append(new_exon)
                        cds_idx += 1
                        break
            ################################################
            # Step 2: parse the inner exons and CDSs
            ################################################
            # Handle the CDSs in the middle
            # Leave the last CDS not processed.
            while cds_idx < len(cds_list)-1:
                cds = cds_list[cds_idx]
                # Create a new exon
                # Add the CDS in the new exon
                new_exon = copy.deepcopy(exon)
                new_exon.update_exon_info(cds.entry.start, cds.entry.end)
                new_exon.add_lifton_cds(cds)
                new_exons.append(new_exon)
                cds_idx += 1
            ################################################
            # Step 3: parse the last CDS & its overlapping exon
            ################################################
            last_cds_processed = False
            cds = cds_list[cds_idx]
            cds_start, cds_end = cds_list[cds_idx].entry.start, cds_list[cds_idx].entry.end
            while exon_idx < len(self.exons):
                exon = self.exons[exon_idx]
                exon_start, exon_end = exon.entry.start, exon.entry.end
                _, ovp = lifton_utils.segments_overlap_length((exon_start, exon_end), (cds_start, cds_end))
                # |eeeeee| |ccccc|
                if exon_end < cds_start:
                    # exon_idx = self.mv_exon_idx(exon_idx)
                    exon_idx += 1
                    continue
                elif ovp:
                    # Create a new exon using exon boundary
                    last_cds_processed = True
                    new_exon = copy.deepcopy(exon)
                    new_exon.update_exon_info(cds_start, max(cds_end, exon_end))
                    new_exon.add_lifton_cds(cds)
                    new_exons.append(new_exon)
                    exon_idx += 1
                    break
                # |cccc|  |eeee|
                elif exon_start > cds_end: 
                    # Create a new exon using CDS boundary
                    last_cds_processed = True
                    new_exon = copy.deepcopy(exon)
                    new_exon.update_exon_info(cds_start, cds_end)
                    new_exon.add_lifton_cds(cds)
                    new_exons.append(new_exon)
                    break
            ################################################
            # Step 4: parse the remaining exons
            ################################################
            while exon_idx < len(self.exons):
                exon = self.exons[exon_idx]
                # Create a new exon
                # Add the CDS in the new exon
                new_exon = copy.deepcopy(exon)
                new_exon.update_exon_info(new_exon.entry.start, new_exon.entry.end)
                new_exons.append(new_exon)
                exon_idx += 1
            ################################################
            # Step 5: parse the last CDS if it does not overlap with any exon
            ################################################
            if not last_cds_processed:
                # Create a new exon
                # Add the CDS in the new exon
                new_exon = copy.deepcopy(exon)
                new_exon.update_exon_info(cds.entry.start, cds.entry.end)
                new_exon.add_lifton_cds(cds)
                new_exons.append(new_exon)
        self.exons = new_exons
        self.update_boundaries()

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
        lcl_exons = []
        lcl_exons = self.exons
        if len(self.exons) > 0 and self.exons[0].entry.strand == '-':
            lcl_exons = self.exons[::-1]
        for exon in lcl_exons:
            # Chaining the exon features
            p_trans_seq = exon.entry.sequence(fai)
            p_trans_seq = Seq(p_trans_seq).upper()
            trans_seq = trans_seq + p_trans_seq
            if exon.cds is not None:
                # Updating CDS frame
                exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)
                # Chaining the CDS features
                p_seq = exon.cds.entry.sequence(fai)
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
        if ref_protein_seq == "" or ref_protein_seq == None:
            lifton_aa_aln = None
            peps = None
        elif protein_seq == "" or protein_seq == None:
            lifton_aa_aln = None
            peps = None
        else:
            peps = protein_seq.split("*")
            lifton_aa_aln = align.protein_align(protein_seq, ref_protein_seq)
            # Update lifton sequence identity
            lifton_status.lifton_aa = max(lifton_status.lifton_aa, lifton_aa_aln.identity)
        return lifton_aa_aln, peps

    def align_trans_seq(self, trans_seq, ref_trans_seq, lifton_status):
        if ref_trans_seq == "" or ref_trans_seq == None:
            lifton_tran_aln = None
        elif trans_seq == "" or trans_seq == None:
            lifton_tran_aln = None
        else:
            lifton_tran_aln = align.trans_align(trans_seq, ref_trans_seq)
            lifton_status.lifton_dna = lifton_tran_aln.identity
        return lifton_tran_aln

    def orf_search_protein(self, fai, ref_protein_seq, ref_trans_seq, lifton_status, is_non_coding, eval_only=False):
        coding_seq, trans_seq = self.get_coding_trans_seq(fai)
        protein_seq = self.translate_coding_seq(coding_seq)
        # Aligning the LiftOn protein & DNA sequences
        lifton_aa_aln, peps = self.align_coding_seq(protein_seq, ref_protein_seq, lifton_status)
        lifton_tran_aln = self.align_trans_seq(trans_seq, ref_trans_seq, lifton_status)
        variants.find_variants(lifton_tran_aln, lifton_aa_aln, lifton_status, peps, is_non_coding)
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
        if ORF_search and eval_only==False:
            self.__find_orfs(trans_seq, ref_protein_seq, lifton_aa_aln, lifton_status)
        return lifton_tran_aln, lifton_aa_aln

    def __find_orfs(self, trans_seq, ref_protein_seq, lifton_aln, lifton_status):
        trans_seq = trans_seq.upper()
        # Find ORFs in whole transcript region (including UTRs)
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
        final_orf = None
        update_orf = False
        max_identity = 0
        for i, orf in enumerate(orf_list):
            orf_DNA_seq = trans_seq[orf.start:orf.end]
            orf_protein_seq = str(Seq(orf_DNA_seq).translate())
            orf_parasail_res = align.parasail_align_protein_base(orf_protein_seq, ref_protein_seq)
            orf_matches, orf_length = get_id_fraction.get_AA_id_fraction(orf_parasail_res.traceback.ref, orf_parasail_res.traceback.query)
            orf_identity = orf_matches/orf_length
            if orf_identity > max_identity:
                max_identity = orf_identity
                final_orf = orf
        # Only update orf if at least one frame similarity is larger than the original lifton_aa by the threshold
        threshold_orf = 0.01
        if max_identity > (lifton_status.lifton_aa + threshold_orf):
            lifton_status.lifton_aa = max_identity
            update_orf = True
        if final_orf is not None and update_orf:
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
            if accum_exon_length <= final_orf.start:
                if final_orf.start < accum_exon_length+curr_exon_len:
                    # Create first partial CDS
                    if exon.cds is not None:
                        if strand == "+":
                            exon.cds.entry.end = exon.entry.end
                            exon.cds.entry.start = exon.entry.start + (final_orf.start - accum_exon_length)
                        elif strand == "-":
                            exon.cds.entry.start = exon.entry.start 
                            exon.cds.entry.end = exon.entry.end - (final_orf.start - accum_exon_length)
                    else:
                        if strand == "+":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.start + (final_orf.start - accum_exon_length), exon.entry.end)
                        elif strand == "-":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.start, exon.entry.end - (final_orf.start - accum_exon_length))
                    exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                    accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)
                else:
                    # No CDS should be created
                    exon.cds = None
            elif final_orf.start < accum_exon_length and accum_exon_length < final_orf.end:
                if final_orf.end <= accum_exon_length+curr_exon_len:
                    # Create the last partial CDS
                    if exon.cds is not None:
                        if strand == "+":
                            exon.cds.entry.start = exon.entry.start
                            exon.cds.entry.end = exon.entry.start + (final_orf.end - accum_exon_length)-1
                        elif strand == "-":
                            exon.cds.entry.end = exon.entry.end
                            exon.cds.entry.start = exon.entry.end - (final_orf.end - accum_exon_length)+1
                    else:
                        if strand == "+":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.start, exon.entry.start + (final_orf.end - accum_exon_length)-1)
                        elif strand == "-":
                            exon.add_novel_lifton_cds(exon.entry, exon.entry.end - (final_orf.end - accum_exon_length)+1, exon.entry.end)
                else:
                    # Keep the original full CDS / extend the CDS to full exon length
                    if exon.cds is None:
                        exon.add_novel_lifton_cds(exon.entry, exon.entry.start, exon.entry.end)
                    else:
                        exon.cds.update_CDS_info(exon.entry.start, exon.entry.end)
                exon.cds.entry.frame = str(self.__get_cds_frame(accum_cds_length))
                accum_cds_length += (exon.cds.entry.end - exon.cds.entry.start + 1)
            elif final_orf.end <= accum_exon_length:
                # No CDS should be created
                exon.cds = None
            accum_exon_length += curr_exon_len

    def __get_cds_frame(self, accum_cds_length):
        return (3 - accum_cds_length%3)%3

    def write_entry(self, fw):
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

    def update_exon_info(self, start, end):
        self.cds = None
        self.entry.source = "LiftOn"
        self.entry.start = start
        self.entry.end = end

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

    def update_CDS_info(self, start, end):
        self.entry.source = "LiftOn"
        self.entry.start = start
        self.entry.end = end

    def write_entry(self, fw):
        fw.write(str(self.entry) + "\n")

    def print_cds(self):
        print(f"\t\t{self.entry}")