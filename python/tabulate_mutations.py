#!/usr/bin/env python

import argparse
from collections import OrderedDict
import os

from Bio import SeqIO


class MutRecord:
    def __init__(self, id, muts):
        self.id = id
        self.muts = muts

    def __repr__(self):
        return "{},{}".format(self.id, " ".join(self.muts))


def find_muts(orig, mutated):
     return [
        "{}{}{}".format(o, idx+1, m)
        for idx, (o, m) in enumerate(zip(orig, mutated)) if o != m]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate sequence differences relative to a naive sequence.")
    parser.add_argument(
        'fasta_path', type=str,
        help="Path to FASTA input alignment.")
    parser.add_argument(
        '--naive', default="naive0",
        help="The name of the naive sequence.")

    args = parser.parse_args()

    seqs = OrderedDict(
        (r.id, str(r.seq)) for r in
        SeqIO.parse(args.fasta_path, "fasta"))

    naive = seqs[args.naive]
    muts = [MutRecord(k, find_muts(naive, s)) for k, s in seqs.items()]

    base, _ = os.path.splitext(args.fasta_path)

    with open(base+'.muts.csv', 'w') as f:
        for mut in muts:
            f.write("{}\n".format(mut))
