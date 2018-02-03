#!/usr/bin/env python

import itertools as it
import argparse
from Bio import SeqIO, Seq, SeqRecord

import ete3
#from ete3 import TreeNode, TreeStyle, NodeStyle, SVG_COLORS
import pickle

ete3 #lint


def load_tree(filename):
    with open(filename, 'rb') as fh:
        return pickle.load(fh)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('simulated_tree_pickle')
    parser.add_argument('sequence_file')
    parser.add_argument('seedid_file')
    return parser.parse_args()


def main():
    args = get_args()
    tree = load_tree(args.simulated_tree_pickle)
    sequences = [SeqRecord.SeqRecord(Seq.Seq(n.sequence), n.name, n.name)
                 for n in it.chain([tree], tree.iter_leaves())]
    seed, _ = tree.get_farthest_leaf()
    with open(args.sequence_file, 'w') as fh:
        SeqIO.write(sequences, fh, 'fasta')
    with open(args.seedid_file, 'w') as fh:
        fh.write(seed.name)

if __name__ == '__main__':
    main()


