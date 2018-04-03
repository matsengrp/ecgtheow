#!/usr/bin/env python

import argparse
import jinja2
import os
import re

from util_functions import parse_fasta_seqs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create a BEAST XML input file from a template.")
    parser.add_argument(
        'template_path', type=str,
        help="Path to XML template.")
    parser.add_argument(
        'fasta_path', type=str,
        help="Path to FASTA input alignment.")
    parser.add_argument(
        '--naive', type=str, required=True,
        help="The name of the naive sequence.")
    parser.add_argument(
        '--seed', type=str, required=True,
        help="The name of the seed sequence.")
    parser.add_argument(
        '--iter', type=int, required=True,
        help="The number of total MCMC iterations.")
    parser.add_argument(
        '--thin', type=int, required=True,
        help="The MCMC sampling frequency.")
    parser.add_argument(
        '--output-dir', type=str,
        help="The name of the output directory.")
    parser.add_argument(
        '--xml-path',
        help="The XML output file path.")

    args = parser.parse_args()

    args.fasta_path = args.fasta_path.lstrip("./")
    if args.xml_path is not None:
        xml_base = os.path.splitext(args.xml_path)[0]
    else:
        fasta_path = re.split("/", args.fasta_path)[-1]
        xml_base = os.path.splitext(fasta_path)[0]
        xml_base = args.output_dir + "/runs/" + xml_base + "_beast"

    id_seq = parse_fasta_seqs(args.fasta_path)

    assert args.naive in id_seq, "Sequence %r not found in FASTA file." % args.naive
    assert args.seed in id_seq, "Sequence %r not found in FASTA file." % args.seed

    temp_vars = dict(
        id_seq=id_seq,
        naive=args.naive,
        iter=args.iter,
        thin=args.thin,
        basename=xml_base
    )

    env = jinja2.Environment(loader = jinja2.FileSystemLoader('.'),
                             undefined=jinja2.StrictUndefined,
                             trim_blocks=True, lstrip_blocks=True)

    env.get_template(args.template_path).stream(**temp_vars).dump(xml_base + ".xml")
