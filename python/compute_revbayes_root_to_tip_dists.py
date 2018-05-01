#!/usr/bin/env python

import argparse
import collections
import csv
import dendropy


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Map trees to root-to-tip distance distributions.")
    parser.add_argument(
        'trees_path', type=str,
        help="Path to RevBayes trees (in BEAST format).")
    parser.add_argument(
        '--burnin', type=int, required=True,
        help="How many entries to remove as burnin.")
    parser.add_argument(
        '--naive', type=str, required=True,
        help="The name of the naive sequence.")
    parser.add_argument(
        '--output-path', type=str, required=True,
        help="The CSV output file path.")

    args = parser.parse_args()

    tree_yielder = dendropy.Tree.yield_from_files(
        files=[args.trees_path],
        schema='nexus',
        preserve_underscores=True
    )
    list_of_dist_lists = []

    for tree_idx, tree in enumerate(tree_yielder):
        if tree_idx < args.burnin:
            # Skip burnin.
            continue

        taxa_labels = [x.label for x in tree.taxon_namespace]
        assert args.naive in taxa_labels, "Please specify the correct naive sequence name."

        root_to_tip_dists = dict()
        tree.calc_node_root_distances(return_leaf_distances_only=True)

        for node in tree.leaf_node_iter():
            if node.taxon.label == args.naive:
                continue
            root_to_tip_dists[node.taxon.label] = node.root_distance + tree.find_node_with_taxon_label(args.naive).root_distance

        root_to_tip_dists = collections.OrderedDict(sorted(root_to_tip_dists.iteritems()))
        list_of_dist_lists.append(root_to_tip_dists)

    keys = [x.keys() for x in list_of_dist_lists]
    assert all(x == keys[0] for x in keys), "Tip names need to be ordered the same."
    values = [x.values() for x in list_of_dist_lists]
    values.insert(0, keys[0])

    with open(args.output_path, "w") as f:
        wr = csv.writer(f)
        wr.writerows(values)
