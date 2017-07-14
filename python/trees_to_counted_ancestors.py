#!/usr/bin/env python

import argparse
from collections import Counter
import math
import numpy as np
from collections import OrderedDict
import os
import sys

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import dendropy
import graphviz

from tabulate_mutations import find_muts


def format_label(label):
    '''
    Format the graphviz edge label to appear "rectangular".
    '''
    label_space = int(np.sqrt(len(label)))
    for i in np.arange(label_space, len(label), label_space + 1):
        label.insert(i, "\n")

    return label


def seqs_of_tree(t, seed):
    '''
    Iterate up the tree, getting ancestral sequences.
    '''
    lineage = [t.find_node_with_taxon_label(seed)]

    while(True):
        n = lineage[-1].parent_node
        if n is None:
            break  # We are done.
        lineage.append(n)

    return [translate(n.annotations.get_value('ancestral')) for n in lineage]


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
    parser.add_argument(
        '--filter', required=True, type=int,
        help="Only display edges with at least this many samples.")

    args = parser.parse_args()
    burnin = int(args.burnin)

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

    source_files = [args.tree_path]
    tree_yielder = dendropy.Tree.yield_from_files(
            files=source_files,
            schema='nexus',
            )

    node_c = Counter()
    edge_c = Counter()
    for tree_idx, t in enumerate(tree_yielder):
        if tree_idx < burnin:
            # Skip burnin.
            continue
        l = seqs_of_tree(t, args.seed)
        l.append(translate(leaf_seqs['naive0']))
        node_c.update(frozenset(l))
        # Note flipped below because we want to go from ancestor to descendant.
        edge_c.update((w, v) for v, w in zip(l[:-1], l[1:]))

    # Count the number of BEAST trees.
    num_trees = node_c.most_common(1)[0][1]

    if node_c.most_common() == []:
        raise Exception ("Nothing to count! Is your burnin too large?")
    out_seqs = OrderedDict()

    # Iterate through all, in order of frequency.
    for idx, (s, count) in enumerate(node_c.most_common(None)):
        # Rename the naive and seed sequences correspondingly if they are in
        # our list.
        if s in special_seqs:
            out_seqs[special_seqs[s]] = s
        else:
            out_seqs['inferred_{}_{}'.format(idx,count)] = s

    if 'naive0' not in out_seqs:
        out_seqs['naive0'] = translate(leaf_seqs['naive0'])

    # Flip the dictionary.
    seqs_out = {v:k for k,v in out_seqs.iteritems()}

    base, _ = os.path.splitext(args.tree_path)

    with open(base+'.aa_lineage_seqs.fasta', 'w') as f:
        for k, v in out_seqs.items():
            f.write('>{}\n'.format(k))
            f.write('{}\n'.format(v))

    dot = graphviz.Digraph(comment=" ".join(sys.argv), format='png')
    dot.attr(size='24,14', ratio='fill', fontsize='14')

# Commented because it defeats our filtering mechanism below. Could re-add if
# we want to add extra information to nodes.
#    for k in out_seqs:
#        dot.node(k)

    for ((a,b), count) in edge_c.most_common(None):
        if a != b and count >= args.filter:
            # Edge confidence measured by percentage of transitions from parent node (i.e. in [0,100]),
            # which is then mapped to the interval [20,100] to avoid transparent edges.
            # Node confidence is treated in a similar fashion below.
            edge_conf = int(20 + (100-20) * float(count) / node_c[a])
            dot.edge(seqs_out[a], seqs_out[b], xlabel=" ".join(format_label(find_muts(a, b))), color="#0000ff" + (str(edge_conf) if edge_conf < 100 else ""), fontsize='11')

            if seqs_out[a] != "naive0":
                child_conf = int(10 + (100-10) * float(node_c[a]) / num_trees)
                dot.node(seqs_out[a], style="filled", fillcolor="#ff0000" + (str(child_conf) if child_conf < 100 else ""))
            if seqs_out[b] != args.seed:
                child_conf = int(10 + (100-10) * float(node_c[b]) / num_trees)
                dot.node(seqs_out[b], style="filled", fillcolor="#ff0000" + (str(child_conf) if child_conf < 100 else ""))

    dot.save(base+'.aa_lineage_graph.dot')
    dot.render(base+'.aa_lineage_graph')
