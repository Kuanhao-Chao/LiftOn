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



class Lifton_GENE:
    def __init__(self, gffutil_entry_gene):
        self.entry = gffutil_entry_gene
        self.transcripts = {}
        
    def add_transcript(self, gffutil_entry_trans):
        # print(">> gffutil_entry_trans[ID]: ", gffutil_entry_trans["ID"])
        Lifton_trans = Lifton_TRANS(gffutil_entry_trans)
        self.transcripts[gffutil_entry_trans["ID"][0]] = Lifton_trans

    def add_exon(self, trans_id, gffutil_entry_exon):
        self.transcripts[trans_id].add_exon(gffutil_entry_exon)

    def add_cds(self, trans_id, gffutil_entry_cds):
        self.transcripts[trans_id].add_cds(gffutil_entry_cds)

    def write_entry(self, fw):
        # print(self.entry)
        fw.write(str(self.entry) + "\n")

    def update_boundaries(self):        
        for key, trans in self.transcripts.items():
            self.entry.start = trans.start if trans.start < self.entry.start else self.entry.start
            self.entry.end = trans.end if trans.end > self.entry.end else self.entry.end


    def print_gene(self):
        print(self.entry)
        for key, trans in self.transcripts.items():
            trans.print_transcript()
        print("\n\n")


class Lifton_TRANS:
    def __init__(self, gffutil_entry_trans):
        self.entry = gffutil_entry_trans
        self.exons = []
        self.exon_dic = {}

    def add_exon(self, gffutil_entry_exon):
        Lifton_exon = Lifton_EXON(gffutil_entry_exon)
        self.exons.append(Lifton_exon)
        self.exon_dic[gffutil_entry_exon.start] = Lifton_exon
        self.exon_dic[gffutil_entry_exon.end] = Lifton_exon

    def add_cds(self, gffutil_entry_cds):
        # self.exons[0].add_cds(gffutil_entry_cds)
        Lifton_exon_retrieval = None
        
        if gffutil_entry_cds.start in self.exon_dic.keys():
            Lifton_exon_retrieval = self.exon_dic[gffutil_entry_cds.start]
        elif gffutil_entry_cds.end in self.exon_dic.keys():
            Lifton_exon_retrieval = self.exon_dic[gffutil_entry_cds.end]

        if Lifton_exon_retrieval is not None:
            Lifton_exon_retrieval.add_cds(gffutil_entry_cds)

    def write_entry(self, fw):
        fw.write(str(self.entry) + "\n")
        
        # Write out the exons first
        for exon in self.exons:
            exon.write_entry(fw)
        # Write out the CDSs first
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
        self.entry = gffutil_entry_exon
        self.cds = None

    def add_cds(self, gffutil_entry_cds):
        Lifton_cds = Lifton_CDS(gffutil_entry_cds)
        self.cds = Lifton_cds

    def write_entry(self, fw):
        # print(self.entry)
        fw.write(str(self.entry) + "\n")

    def print_exon(self):
        print(f"\t\t{self.entry}")
        if self.cds != None:
            self.cds.print_cds()


class Lifton_CDS:
    def __init__(self, gffutil_entry_cds):
        self.entry = gffutil_entry_cds

    def write_entry(self, fw):
        # print(self.entry)
        fw.write(str(self.entry) + "\n")

    def print_cds(self):
        print(f"\t\t{self.entry}")
