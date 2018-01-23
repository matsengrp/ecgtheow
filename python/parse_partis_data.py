#!/usr/bin/env python

import argparse
import os
import yaml


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse partis cluster data and generate 'healthy' seeded sequences.")
    parser.add_argument(
        'data_dir', type=str,
        help="Path to partis data directory.")
    parser.add_argument(
        '--sample', type=str, required=True,
        help="The name of the sample.")
    parser.add_argument(
        '--seed', type=str, required=True,
        help="The name of the seed sequence.")

    args = parser.parse_args()
    args.data_dir = args.data_dir.rstrip("/")

    # Load the YAML file.
    with open(args.data_dir + '/info.yaml', 'r') as f:
        yaml_dict = yaml.load(f)

    param_dir = yaml_dict["samples"][args.sample]["parameter-dir"]
    locus = yaml_dict["samples"][args.sample]["meta"]["locus"]
    annotation_path = yaml_dict["samples"][args.sample]["seeds"][args.seed]["cluster-annotation-file"]
    partition_path = yaml_dict["samples"][args.sample]["seeds"][args.seed]["partition-file"]

    os.system("lib/cft/bin/process_partis.py" +\
              " --cluster-annotation-file " + annotation_path +\
              " --partition-file " + partition_path +\
              " --seqs-out " + "data/" + args.seed + ".family_0.healthy.fasta" +\
              " --parameter-dir " + param_dir +\
              " --locus " + locus +\
              " --remove-frameshifts --remove-stops" +\
              " --remove-mutated-invariants --indel-reversed-seqs")
