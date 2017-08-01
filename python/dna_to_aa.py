#!/usr/bin/env python

import argparse
import glob
from collections import OrderedDict
import os

from util_functions import parse_fasta_seqs, translate, write_to_fasta


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert DNA sequences to amino acid sequences in a FASTA file.")
    parser.add_argument(
        'fasta_path', type=str,
        help="Path to a FASTA file.")

    args = parser.parse_args()

    id_seq = parse_fasta_seqs(args.fasta_path)
    id_seq = OrderedDict([(k,translate(v)) for (k,v) in id_seq.iteritems()])

    base, ext = os.path.splitext(args.fasta_path)
    write_to_fasta(id_seq, base + ".AA" + ext)
