import parasail

class Alignment():
    def __init__(self, ref_seq, target_seq):
        self.ref_seq = ref_seq
        self.target_seq = target_seq
        self.alignment = self.align_seqs()


    def calculate_seq_id(self, matches=None, len=None):
        if matches is None and len  is None:
            matches, len = self.get_id_fraction()
        return str(matches / len)[0:5]


    def get_id_fraction(self):
        matches = 0
        seq1 = self.alignment.traceback.query
        seq2 = self.alignment.traceback.ref
        for i, letter in enumerate(seq1):
            if letter == seq2[i]:
                matches += 1
        return matches, max(len(seq1), len(seq2))



class ProteinAlignment(Alignment):
    def __init__(self, ref_seq, target_seq):
        super().__init__(ref_seq, target_seq)


    def align_seqs(self):
        matrix = parasail.Matrix("blosum62")
        gap_open = 11
        gap_extend = 1
        return parasail.nw_trace_scan_sat( self.ref_seq, self.target_seq, gap_open,gap_extend,matrix)


class DNAAlignment(Alignment):
    def __init__(self, ref_seq, target_seq):
        super().__init__(ref_seq, target_seq)


    def align_seqs(self):
        matrix = parasail.matrix_create("ACGT*", 1, -3)
        gap_open = 5
        gap_extend = 2
        return parasail.nw_trace_scan_sat(self.ref_seq, self.target_seq,gap_open,gap_extend,matrix)
