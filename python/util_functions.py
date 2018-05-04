
from collections import OrderedDict

from Bio.Alphabet import Gapped, IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

def parse_fasta_seqs(path, invert=False):
    '''
    Parse a FASTA file and create a (id:seq) or (seq:id) ordered dict.
    '''
    return OrderedDict(
        (str(r.seq), r.id) if invert else (r.id, str(r.seq))
        for r in SeqIO.parse(path, "fasta")
    )

def seqs_of_tree(t, seed):
    '''
    Iterate up the tree, getting ancestral sequences.
    '''
    lineage = [t.find_node_with_taxon_label(seed)]

    while(True):
        node = lineage[-1].parent_node
        if node is None:
            break  # We are done.
        lineage.append(node)

    return [n.annotations.get_value('ancestral') for n in lineage]

def translate(s):
    '''
    Assume we are in frame and translate DNA to amino acids.
    '''
    coding_dna = Seq(s[:(3*int(len(s)/3))], Gapped(IUPAC.ambiguous_dna))
    return str(coding_dna.translate())

def write_to_fasta(d, filename):
    '''
    Write a FASTA dict to file.
    '''
    with open(filename, 'w') as f:
        for k, v in d.items():
            f.write('>{}\n'.format(k))
            f.write('{}\n'.format(v))
        f.close()
