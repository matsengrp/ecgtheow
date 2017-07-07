#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
import dendropy
import os
import sys


def seqs_of_tree(t, seed):
    '''
    Iterate up the tree, looking at ancestral sequences.
    '''
    lineage = [t.find_node_with_taxon_label(seed)]

    while(True):
        n = lineage[-1].parent_node
        if n is None:
            break
        lineage.append(n)

    return frozenset([
        translate(n.annotations.get_value('ancestral')) for n in lineage])


def translate(s):
    '''
    Assume we are in frame and translate to amino acids
    '''
    coding_dna = Seq(s[:(3*int(len(s)/3))], IUPAC.unambiguous_dna)
    return str(coding_dna.translate())


in_name = sys.argv[1]
burnin = 20
source_files = [in_name]
tree_yielder = dendropy.Tree.yield_from_files(
        files=source_files,
        schema='nexus',
        )

c = Counter()
for tree_idx, t in enumerate(tree_yielder):
    if tree_idx < burnin:
        # Skip burnin.
        continue
    c.update(seqs_of_tree(t, 'BF520.1-igh'))


base, _ = os.path.splitext(in_name)

with open(base+'.aa_lineage.fasta', 'w') as f:
    for idx, (s, count) in enumerate(c.most_common(None)):
        f.write('>{}_{}\n'.format(idx,count))
        f.write(s+'\n')

