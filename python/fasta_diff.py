#!/usr/bin/env python

import argparse
import os
import re

from util_functions import parse_fasta_seqs


def filter_beast_seqs(seqs, seq_id, rgx, nfilter):
    filtered_seqs = set()
    for seq in seqs:
        check = rgx.findall(seq_id[seq])
        if not check or (check and int(check[0]) >= nfilter):
            filtered_seqs.add(seq)

    return filtered_seqs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compare two FASTA files containing BEAST ancestral sequence counts.")
    parser.add_argument(
        'file1_path', type=str,
        help="Path to a FASTA file.")
    parser.add_argument(
        'file2_path', type=str,
        help="Path to a FASTA file.")
    parser.add_argument(
        '--filter', type=int, default=0,
        help="The threshold used to remove infrequent ancestral sequences.")

    args = parser.parse_args()

    seq1_id1 = parse_fasta_seqs(args.file1_path, invert=True)
    seq2_id2 = parse_fasta_seqs(args.file2_path, invert=True)

    rgx = re.compile("^inferred_[0-9]+_([0-9]+)$")

    seq1_keys = set(seq1_id1)
    seq2_keys = set(seq2_id2)
    seq1_keys = filter_beast_seqs(seq1_keys, seq1_id1, rgx, args.filter)
    seq2_keys = filter_beast_seqs(seq2_keys, seq2_id2, rgx, args.filter)

    seq12_set = seq1_keys & seq2_keys
    seq1_set = seq1_keys - seq2_keys
    seq2_set = seq2_keys - seq1_keys

    print "FILE1 & FILE2:\n"

    for elem in seq12_set:
        print seq1_id1[elem] + " <---> " + seq2_id2[elem]

    print "\nFILE1 ONLY:\n"

    for elem in seq1_set:
        print seq1_id1[elem]

    print "\nFILE2 ONLY:\n"

    for elem in seq2_set:
        print seq2_id2[elem]
