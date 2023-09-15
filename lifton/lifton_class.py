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
        Lifton_trans = Lifton_TRANS(gffutil_entry_trans)
        
        print(">> gffutil_entry_trans[ID]: ", gffutil_entry_trans["ID"])
        self.transcripts[gffutil_entry_trans["ID"][0]] = Lifton_trans

    def add_transcript_exon(self, trans_id, gffutil_entry_exon):
        self.transcripts[trans_id].add_exon(gffutil_entry_exon)

    def add_transcript_cds(self, trans_id, gffutil_entry_cds):
        pass

    def write_entry(self, fw):
        print(self.entry)
        # fw.write(str(self.entry) + "\n")

    def update_boundaries(self):
        self.entry.start = 0
        self.entry.end = 1
        # pass
        # print(self.entry)
        # fw.write(str(self.entry) + "\n")

    def print(self):
        print(f'Feature ID: {self.entry.id}')
        print(f'Chromosome: {self.entry.chrom}')
        print(f'Start: {self.entry.start}')
        print(f'End: {self.entry.end}')
        print(f'Strand: {self.entry.strand}')



class Lifton_TRANS:
    def __init__(self, gffutil_entry_trans):
        self.entry = gffutil_entry_trans
        self.exons = []

    def add_exon(self, gffutil_entry_exon):
        Lifton_exon = Lifton_EXON(gffutil_entry_exon)
        self.exons.append(Lifton_exon)

    def add_cds(self, gffutil_entry_cds):
        self.exons[0].add_cds(gffutil_entry_cds)


class Lifton_EXON:
    def __init__(self, gffutil_entry_exon):
        self.entry = gffutil_entry_exon
        self.cds = None

    def add_cds(self, gffutil_entry_cds):
        Lifton_cds = Lifton_CDS(gffutil_entry_cds)
        self.cds = Lifton_cds



class Lifton_CDS:
    def __init__(self, gffutil_entry_cds):
        self.entry = gffutil_entry_cds
