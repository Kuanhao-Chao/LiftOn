from lifton import lifton_utils
import copy

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
        gffutil_entry_gene.source = "Lifton"
        
    def add_transcript(self, gffutil_entry_trans):
        # print(">> gffutil_entry_trans[ID]: ", gffutil_entry_trans["ID"])
        Lifton_trans = Lifton_TRANS(gffutil_entry_trans)
        self.transcripts[gffutil_entry_trans["ID"][0]] = Lifton_trans
        return Lifton_trans

    def add_exon(self, trans_id, gffutil_entry_exon):
        self.transcripts[trans_id].add_exon(gffutil_entry_exon)

    def add_cds(self, trans_id, gffutil_entry_cds):
        self.transcripts[trans_id].add_cds(gffutil_entry_cds)
                            
                            
    def update_cds_list(self, trans_id, cds_list):
        # print("Inside 'self.transcripts[trans_id]': ", self.transcripts[trans_id])
        trans_selected = self.transcripts[trans_id]
        trans_selected.update_cds_list(cds_list)
        self.transcripts[trans_id] = trans_selected
     
        self.update_boundaries()
        # print("\tnew cds_list (inner): ", len(self.transcripts[trans_id].exons))


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
        gffutil_entry_trans.source = "Lifton"

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

    def update_cds_list(self, cds_list):
        print(f"\t>> update_cds_list (len: {len(cds_list)}) ")
        print(f"\t>> self.exons (len: {len(self.exons)}) ")

        idx_exon_itr = 0

        first_cds = cds_list[0]
        last_cds = cds_list[-1]

        new_exons = []

        # print("len(self.exons): ", len(self.exons))

        while idx_exon_itr < len(self.exons)-1:
            exon = self.exons[idx_exon_itr]
            # print("> exon: ", exon)
            # print("> first_cds: ", first_cds)
            # print(f"Checking overlapping: {exon.entry.start}-{exon.entry.end}; {first_cds.entry.start}-{first_cds.entry.end}")
            
            if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (first_cds.entry.start, first_cds.entry.end)):
                # print("## Overlap!!!")
                # 1. Create a new exon
                # 2. Create a new CDS

                # new_exon = copy.deepcopy(exon)
                # new_cds = copy.deepcopy(first_cds)

                if exon.entry.start <= first_cds.entry.start:
                    exon.entry.end = first_cds.entry.end
                elif exon.entry.start > first_cds.entry.start:
                    exon.entry.start = first_cds.entry.start
                    exon.entry.end = first_cds.entry.end
                
                exon.add_lifton_cds(first_cds)
                new_exons.append(exon)

                break
            
            elif exon.entry.end <= first_cds.entry.start:
                exon.add_lifton_cds(None)
                new_exons.append(exon)
            elif exon.entry.start >= last_cds.entry.end:
                break
            idx_exon_itr += 1


        for inner_cds in cds_list[1:len(cds_list)-1]:
            # 1. Create a new exon
            # 2. Create a new CDS
            new_inner_exon = copy.deepcopy(exon)
            new_inner_exon.entry.start = inner_cds.entry.start
            new_inner_exon.entry.end = inner_cds.entry.end
            new_inner_exon.add_lifton_cds(inner_cds)
            new_exons.append(new_inner_exon)


        while idx_exon_itr < len(self.exons):
            # print("idx_exon_itr: ", idx_exon_itr)
            exon = self.exons[idx_exon_itr]
            
            if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (last_cds.entry.start, last_cds.entry.end)):
                # 1. Create a new exon
                # 2. Create a new CDS
                if exon.entry.end <= last_cds.entry.end:
                    exon.entry.start = last_cds.entry.start
                    exon.entry.end = last_cds.entry.end

                elif exon.entry.end > last_cds.entry.end:
                    exon.entry.start = last_cds.entry.start
                
                exon.add_lifton_cds(last_cds)
                new_exons.append(exon)

            elif exon.entry.end <= last_cds.entry.start:
                # just skip. CDSs have been already created.
                pass
            elif exon.entry.start >= last_cds.entry.end:
                # create exon only
                exon.add_lifton_cds(None)
                new_exons.append(exon)


            idx_exon_itr += 1

        print(f"\t>> new_exons (len: {len(new_exons)}) ")

        # for n_exon in new_exons:
        #     n_exon.print_exon()
        self.exons = new_exons


    def write_entry(self, fw):
        # print("Inside 'write_entry'!")
        # print(self.entry)
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
        gffutil_entry_exon.source = "Lifton"

    def add_cds(self, gffutil_entry_cds):
        Lifton_cds = Lifton_CDS(gffutil_entry_cds)
        self.cds = Lifton_cds

    def add_lifton_cds(self, Lifton_cds):
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
        gffutil_entry_cds.source = "Lifton"

    def write_entry(self, fw):
        # print(self.entry)
        fw.write(str(self.entry) + "\n")

    def print_cds(self):
        print(f"\t\t{self.entry}")
