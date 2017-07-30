
from collections import OrderedDict
from Bio import SeqIO

def parse_fasta_seqs(path, invert=False):
    '''
    Parse a FASTA file and create a (id:seq) or (seq:id) ordered dict.
    '''
    return OrderedDict(
        (str(r.seq), r.id) if invert else (r.id, str(r.seq))
        for r in SeqIO.parse(path, "fasta")
    )
