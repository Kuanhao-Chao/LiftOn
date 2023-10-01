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
        gffutil_entry_gene.source = "Lifton"
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
                            
                            
    def update_cds_list(self, trans_id, cds_list):
        self.transcripts[trans_id].update_cds_list(cds_list)
        self.update_boundaries()
        # print("\tnew cds_list (inner): ", len(self.transcripts[trans_id].exons))     


    def write_entry(self, fw):
        # print(self.entry)
        fw.write(str(self.entry) + "\n")
        for key, trans in self.transcripts.items():
            trans.write_entry(fw)

    def update_boundaries(self):        
        # print("\tself.transcripts length: ", len(self.transcripts))
        for key, trans in self.transcripts.items():
            # print("\t## key: ", key)
            self.entry.start = trans.entry.start if trans.entry.start < self.entry.start else self.entry.start
            self.entry.end = trans.entry.end if trans.entry.end > self.entry.end else self.entry.end

        # print(f"update_boundaries:  {self.entry.start}-{self.entry.end}")

    def print_gene(self):
        # print(self.entry)
        for key, trans in self.transcripts.items():
            trans.print_transcript()
        # print("\n\n")


class Lifton_TRANS:
    def __init__(self, gffutil_entry_trans):
        gffutil_entry_trans.source = "Lifton"
        self.entry = gffutil_entry_trans
        self.exons = []
        self.exon_dic = {}

    def add_exon(self, gffutil_entry_exon):
        # print("!!@>> Run add_exon function! ")
        Lifton_exon = Lifton_EXON(gffutil_entry_exon)
        lifton_utils.custom_bisect_insert(self.exons, Lifton_exon)
        # self.exons.append(Lifton_exon)
        # self.exon_dic[gffutil_entry_exon.start] = Lifton_exon
        # self.exon_dic[gffutil_entry_exon.end] = Lifton_exon

    def add_cds(self, gffutil_entry_cds):
        # print("!!@>> Run add_cds function! ")
        # # self.exons[0].add_cds(gffutil_entry_cds)
        # Lifton_exon_retrieval = None
        # if gffutil_entry_cds.start in self.exon_dic.keys():
        #     Lifton_exon_retrieval = self.exon_dic[gffutil_entry_cds.start]
        # elif gffutil_entry_cds.end in self.exon_dic.keys():
        #     Lifton_exon_retrieval = self.exon_dic[gffutil_entry_cds.end]

        for exon in self.exons:
            if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (gffutil_entry_cds.start, gffutil_entry_cds.end)):
                exon.add_cds(gffutil_entry_cds)

    def update_gffutil_entry_trans(self, gffutil_entry_trans):
        print("gffutil_entry_trans: ", gffutil_entry_trans.attributes)

        for key, atr in gffutil_entry_trans.attributes.items():
            self.entry.attributes[key] = atr

    def update_cds_list(self, cds_list):
        print(f"\t>> update_cds_list (len: {len(cds_list)}) ")
        print(f"\t>> self.exons (len: {len(self.exons)}) ")

        idx_exon_itr = 0
        new_exons = []

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
                # print("> exon: ", exon)
                # print("> first_cds: ", first_cds)
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
            
            # while idx_exon_itr < len(self.exons):
            #     exon = self.exons[idx_exon_itr]
            #     idx_exon_itr += 1

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

            # cds.print_cds()
            # exon.print_exon()

            # last_exon = self.exons[-1]
            # last_cds = cds_list[-1]

            # print("last_cds: ")
            # last_cds.print_cds()

            # print(f"All cdss ({len(cds_list)}): ")
            # for cds_p in cds_list:
            #     cds_p.print_cds()


            # print(f"All exons ({len(self.exons)}): ")
            # for exon_p in self.exons:
            #     exon_p.print_exon()



            # while cds_idx < len(cds_list) or exon_idx < len(self.exons):
            #     if cds_idx == len(cds_list):
            #         exon_idx += 1
            #     elif exon_idx == len(self.exons):
            #         cds_idx += 1
                
            #     exon = self.exons[exon_idx]
            #     cds = cds_list[cds_idx]

            #     if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (cds.entry.start, cds.entry.end)):
            #         pass
            #     else:

            ################################################
            # Step 1: CDSs or exons are smaller & no overlapping 
            #      => finding the first overlapping CDS and exons
            #      => processing all exons / CDS ahead of the first overlapping
            ################################################
            # print("\n\n")

            init_head_order = None
            if exon.entry.start >= cds.entry.end:
                init_head_order = 0
            elif exon.entry.end <= cds.entry.start:
                init_head_order = 1
                
            while not lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (cds.entry.start, cds.entry.end)):
                print(f"cds_idx: {cds_idx}; {cds.entry.start}-{cds.entry.end} (len: {len(cds_list)});  exon_idx: {exon_idx}; {exon.entry.start}-{exon.entry.end} (len: {len(self.exons)})")
                if exon.entry.start >= cds.entry.end:

                    # print("Order: ", 0)
                    if init_head_order != 0:
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
            # print(">> parse the first overlapping exon & CDS")
            # print(f"cds_idx: {cds_idx}; {cds.entry.start}-{cds.entry.end} (len: {len(cds_list)});  exon_idx: {exon_idx}; {exon.entry.start}-{exon.entry.end} (len: {len(self.exons)})")

            new_exon = copy.deepcopy(exon)
            if exon.entry.start > cds.entry.start:
                new_exon.entry.start = cds.entry.start
            if exon.entry.end is not cds.entry.end:
                new_exon.entry.end = cds.entry.end
            
            new_exon.add_lifton_cds(cds)
            new_exons.append(new_exon)
            cds_idx += 1

            ################################################
            # Step 3: parse the inner exons and CDSs
            ################################################
            # Handle the CDSs in the middle
            # Leave the last CDS not processed.
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

            # print(">>>>>>>>>>>START!!!")
            # print("exon_idx: ", exon_idx)
            # exon.print_exon()
            # print("cds_idx: ", cds_idx)
            # cds.print_cds()
            
            while exon_idx < len(self.exons):
                exon = self.exons[exon_idx]
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
                        print(">> !!New exon!!")
                        new_exon.add_lifton_cds(cds)
                        new_exons.append(new_exon)
                # new_exon.print_exon()
                exon_idx += 1


            # ################################################
            # # Step 5: parse all remaining exons without CDS if there's any
            # ################################################


            # # # print("len(self.exons): ", len(self.exons))
            # # # Handle the first CDS            
            # # while idx_exon_itr < len(self.exons):
            # #     exon = self.exons[idx_exon_itr]
            # #     idx_exon_itr += 1
            # #     # print("> exon: ", exon)
            # #     # print("> first_cds: ", first_cds)
            # #     # print(f"Checking overlapping: {exon.entry.start}-{exon.entry.end}; {first_cds.entry.start}-{first_cds.entry.end}")
                
            # #     if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (first_cds.entry.start, first_cds.entry.end)):
            # #         # print("## Overlap!!!")
            # #         # 1. Create a new exon
            # #         # 2. Create a new CDS

            # #         # new_exon = copy.deepcopy(exon)
            # #         # new_cds = copy.deepcopy(first_cds)

            # #         if exon.entry.start <= first_cds.entry.start:
            # #             exon.entry.end = first_cds.entry.end
            # #         elif exon.entry.start > first_cds.entry.start:
            # #             exon.entry.start = first_cds.entry.start
            # #             exon.entry.end = first_cds.entry.end
                    
            # #         exon.add_lifton_cds(first_cds)
            # #         new_exons.append(exon)

            # #         break
                
            # #     elif exon.entry.end <= first_cds.entry.start:
            # #         exon.add_lifton_cds(None)
            # #         new_exons.append(exon)
            # #     elif exon.entry.start >= last_cds.entry.end:
            # #         break




            # # Handle the last CDSs
            # # print("idx_exon_itr: ", idx_exon_itr)
            # while idx_exon_itr < len(self.exons):
            #     # print("idx_exon_itr: ", idx_exon_itr)
            #     exon = self.exons[idx_exon_itr]
            #     # exon.print_exon()
                
            #     if lifton_utils.segments_overlap((exon.entry.start, exon.entry.end), (last_cds.entry.start, last_cds.entry.end)):
            #         # print("\tCase 1")
            #         # 1. Create a new exon
            #         # 2. Create a new CDS
            #         if exon.entry.end <= last_cds.entry.end:
            #             exon.entry.start = last_cds.entry.start
            #             exon.entry.end = last_cds.entry.end

            #         elif exon.entry.end > last_cds.entry.end:
            #             exon.entry.start = last_cds.entry.start
                    
            #         exon.add_lifton_cds(last_cds)
            #         new_exons.append(exon)

            #     elif exon.entry.end <= last_cds.entry.start:
            #         # just skip. CDSs have been already created.
            #         # print("\tCase 2")
            #         pass
            #     elif exon.entry.start >= last_cds.entry.end:
            #         # create exon only
            #         # print("\tCase 3")
            #         exon.add_lifton_cds(None)
            #         new_exons.append(exon)


            #     idx_exon_itr += 1

        # print(f"\t>> new_exons (len: {len(new_exons)}) ")


        # print(f"All exons ({len(new_exons)}): ")
        # for exon in new_exons:
        #     exon.print_exon()


        # for n_exon in new_exons:
        #     n_exon.print_exon()
        self.exons = new_exons
        self.update_boundaries()


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
        gffutil_entry_exon.source = "Lifton"
        gffutil_entry_exon.featuretype = "exon"
        self.entry = gffutil_entry_exon
        self.cds = None

        # print("self.entry: ", self.entry.featuretype)

    def update_exon_info(self, start, end):
        self.cds = None
        self.entry.source = "Lifton"
        self.entry.start = start
        self.entry.end = end
        # print(f"start: {self.entry.start} - end: {self.entry.end}")

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
        gffutil_entry_cds.source = "Lifton"
        gffutil_entry_cds.featuretype = "CDS"
        self.entry = gffutil_entry_cds

    def write_entry(self, fw):
        # print(self.entry)
        fw.write(str(self.entry) + "\n")

    def print_cds(self):
        print(f"\t\t{self.entry}")
