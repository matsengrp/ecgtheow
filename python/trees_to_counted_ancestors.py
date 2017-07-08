#!/usr/bin/env python

import argparse
from collections import Counter
from collections import OrderedDict
import os
import sys

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import dendropy


def seqs_of_tree(t, seed):
    '''
    Iterate up the tree, looking at ancestral sequences.
    '''
    lineage = [t.find_node_with_taxon_label(seed)]

    while(True):
        n = lineage[-1].parent_node
        if n is None:
            break  # We are done.
        lineage.append(n)

    return frozenset([
        translate(n.annotations.get_value('ancestral')) for n in lineage])


def translate(s):
    '''
    Assume we are in frame and translate to amino acids
    '''
    coding_dna = Seq(s[:(3*int(len(s)/3))], IUPAC.unambiguous_dna)
    return str(coding_dna.translate())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count occurrences of BEAST-inferred ancestral sequences.")
    parser.add_argument(
        'tree_path', type=str,
        help="Path to BEAST output.")
    parser.add_argument(
        'fasta_path', type=str,
        help="Path to FASTA input alignment.")
    parser.add_argument(
        '--burnin', required=True,
        help="How many entries to remove as burnin.")
    parser.add_argument(
        '--seed', required=True,
        help="The name of the seed sequence.")

    args = parser.parse_args()
    burnin = int(args.burnin)

    source_files = [args.tree_path]
    tree_yielder = dendropy.Tree.yield_from_files(
            files=source_files,
            schema='nexus',
            )

    c = Counter()
    for tree_idx, t in enumerate(tree_yielder):
        if tree_idx < burnin:
            # Skip burnin.
            continue
        c.update(seqs_of_tree(t, args.seed))

    leaf_seqs = {
        k:str(v.seq) for k,v in
        SeqIO.to_dict(SeqIO.parse(args.fasta_path, "fasta")).items() }

    # Make a reversed dictionary for two special sequences.
    special_seqs = {}
    try:
        if leaf_seqs['naive0'] == leaf_seqs[args.seed]:
            raise Exception("The naive sequence is the same as the seed sequence!")
        for special_name in ['naive0', args.seed]:
            special_seqs[translate(leaf_seqs[special_name])] = special_name
    except KeyError, e:
        raise Exception("Couldn't find {} in FASTA file {}.".format(e, args.fasta_path))

    out_seqs = OrderedDict()

    # Iterate through all, in order of frequency.
    for idx, (s, count) in enumerate(c.most_common(None)):
        # Rename the naive and seed sequences correspondingly if they are in our list.
        if s in special_seqs:
            out_seqs[special_seqs[s]] = s
        else:
            out_seqs['inferred_{}_{}'.format(idx,count)] = s

    if 'naive0' not in out_seqs:
        out_seqs['naive0'] = translate(leaf_seqs['naive0'])

    base, _ = os.path.splitext(args.tree_path)

    with open(base+'.aa_lineage.fasta', 'w') as f:
        for k, v in out_seqs.items():
            f.write('>{}\n'.format(k))
            f.write('{}\n'.format(v))
