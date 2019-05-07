#!/usr/bin/env python

import argparse
import os
import yaml

from util_functions import getsuffix


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
    parser.add_argument(
        '--other-partition-id', type=str,
        help="The key of the other partition to use from the seed's other-partitions if they exist.")
    parser.add_argument(
        '--output-path', type=str, required=True,
        help="Path to output FASTA file.")

    args = parser.parse_args()
    args.data_dir = args.data_dir.rstrip("/")

    # Load the YAML file.
    with open(args.data_dir + '/info.yaml', 'r') as f:
        yaml_dict = yaml.load(f)

    sample = yaml_dict["samples"][args.sample]
    locus = sample["meta"]["locus"]
    if args.other_partition_id is not None:
        partition_path = sample["seeds"][args.seed].get("other-partitions", {}).get(args.other_partition_id, {}).get("partition-file")
        if partition_path is None:
            raise Exception('key {} not found in other-partitions for sample: {}, seed: {}. Either because there is no other-partitions key for this seed or because the key specified is not in the other-partitions list. Check info.yaml to confirm you are using the correct string for --other-partition-id for your sample and seed'.format(args.other_partition_id, args.sample, args.seed))
    else:
        partition_path = sample["seeds"][args.seed]["partition-file"]

    os.system("export PARTIS=${PWD%/}/lib/cft/partis;" +\
              "lib/cft/bin/process_partis.py" +\
              " --partition-file " + partition_path +\
              " --seqs-out " + args.output_path +\
             (" --glfo-dir " + sample["glfo-dir"] if getsuffix(partition_path) == ".csv" else "") +\
              " --locus " + locus +\
              " --inferred-naive-name naive" +\
              " --remove-frameshifts --remove-stops" +\
              " --remove-mutated-invariants --indel-reversed-seqs")
