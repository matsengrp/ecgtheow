#!/usr/bin/env python

import argparse
import os
import re
from Bio import SeqIO


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Compare two FASTA files containing BEAST ancestral sequence counts.")
    parser.add_argument(
        'file1_path', type=str,
        help="Path to a FASTA file.")
    parser.add_argument(
        'file2_path', type=str,
        help="Path to a FASTA file.")
    parser.add_argument(
        '--naive', type=str, required=True,
        help="The name of the naive sequence.")
    parser.add_argument(
        '--seed', type=str, required=True,
        help="The name of the seed sequence.")
    parser.add_argument(
        '--filter', type=int, default=0,
        help="The threshold used to remove infrequent ancestral sequences.")

    args = parser.parse_args()

    assert os.path.splitext(args.file1_path)[1] == ".fasta", "%r is not a proper FASTA file path" % args.file1_path
    assert os.path.splitext(args.file2_path)[1] == ".fasta", "%r is not a proper FASTA file path" % args.file2_path

    seq1_id1 = {
        str(v.seq):k for k,v in
        SeqIO.to_dict(SeqIO.parse(args.file1_path, "fasta")).items() }
    id1_seq1 = { v:k for k,v in seq1_id1.iteritems() }
    seq2_id2 = {
        str(v.seq):k for k,v in
        SeqIO.to_dict(SeqIO.parse(args.file2_path, "fasta")).items() }
    id2_seq2 = { v:k for k,v in seq2_id2.iteritems() }

    assert args.naive in id1_seq1 and args.naive in id2_seq2, "Sequence %r not found in both FASTA files." % args.naive
    assert args.seed in id1_seq1 and args.seed in id2_seq2, "Sequence %r not found in both FASTA files." % args.seed

    seq1_keys = set(seq1_id1)
    seq2_keys = set(seq2_id2)
    seq12_set = seq1_keys & seq2_keys
    seq1_set = seq1_keys - seq2_keys
    seq2_set = seq2_keys - seq1_keys

    print "FILE1 && FILE2:\n"

    for elem in seq12_set:
        if re.compile("^inferred_[0-9]+_([0-9]+)$").findall(seq1_id1[elem]):
            if int(re.compile("^inferred_[0-9]+_([0-9]+)$").findall(seq1_id1[elem])[0]) < args.filter:
                continue
        if re.compile("^inferred_[0-9]+_([0-9]+)$").findall(seq2_id2[elem]):
            if int(re.compile("^inferred_[0-9]+_([0-9]+)$").findall(seq2_id2[elem])[0]) < args.filter:
                continue
        if re.compile("^([0-9]+)$").findall(seq1_id1[elem]):
            if int(re.compile("^([0-9]+)$").findall(seq1_id1[elem])[0]) < args.filter:
                continue
        if re.compile("^([0-9]+)$").findall(seq2_id2[elem]):
            if int(re.compile("^([0-9]+)$").findall(seq2_id2[elem])[0]) < args.filter:
                continue
        if seq1_id1[elem] == args.naive or seq1_id1[elem] == args.seed:
            continue

        print seq1_id1[elem] + " <---> " + seq2_id2[elem]

    print "\nFILE1 ONLY:\n"

    for elem in seq1_set:
        if re.compile("^inferred_[0-9]+_([0-9]+)$").findall(seq1_id1[elem]):
            if int(re.compile("^inferred_[0-9]+_([0-9]+)$").findall(seq1_id1[elem])[0]) < args.filter:
                continue
        if re.compile("^([0-9]+)$").findall(seq1_id1[elem]):
            if int(re.compile("^([0-9]+)$").findall(seq1_id1[elem])[0][0]) < args.filter:
                continue

        print seq1_id1[elem]

    print "\nFILE2 ONLY:\n"

    for elem in seq2_set:
        if re.compile("^inferred_[0-9]+_([0-9]+)$").findall(seq2_id2[elem]):
            if int(re.compile("^inferred_[0-9]+_([0-9]+)$").findall(seq2_id2[elem])[0]) < args.filter:
                continue
        if re.compile("^([0-9]+)$").findall(seq2_id2[elem]):
            if int(re.compile("^([0-9]+)$").findall(seq2_id2[elem])[0][0]) < args.filter:
                continue

        print seq2_id2[elem]
