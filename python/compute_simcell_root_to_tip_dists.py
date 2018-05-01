#!/usr/bin/env python

import argparse
import collections
import csv
import ete3
import pickle


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Map simulated ETE tree to root-to-tip distance distribution.")
    parser.add_argument(
        'tree_path', type=str,
        help="Path to simulated ETE tree (as a pickle).")
    parser.add_argument(
        '--output-path', type=str, required=True,
        help="The CSV output file path.")

    args = parser.parse_args()

    with open(args.tree_path) as fh:
        tree = pickle.load(fh)

    list_of_dists = []

    for leaf_node in tree.get_leaves():
        list_of_dists.append(leaf_node.get_distance(tree))

    with open(args.output_path, "w") as f:
        wr = csv.writer(f)
        wr.writerows([list_of_dists])
