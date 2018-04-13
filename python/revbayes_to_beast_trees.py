#!/usr/bin/env python

import argparse
import dendropy
import os
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert RevBayes-annotated trees to BEAST-annotated trees.")
    parser.add_argument(
        'trees_path', type=str,
        help="Path to RevBayes trees file.")
    parser.add_argument(
        'asrs_path', type=str,
        help="Path to RevBayes ASR file.")
    parser.add_argument(
        '--output-path', type=str,
        help="Path to output BEAST trees.")

    args = parser.parse_args()

    output_base = os.path.splitext(args.output_path or args.trees_path)[0]

    # Annotate the trees with ASR sequences.
    tree_strs = pd.read_csv(args.trees_path, sep="\t")["psi"].tolist()
    asrs = pd.read_csv(args.asrs_path, sep="\t", comment="#")

    trees = dendropy.TreeList()
    for i, tree_str in enumerate(tree_strs):
        tree = dendropy.Tree.get(data=tree_str, schema="newick")

        for node in tree.postorder_node_iter():
            node_ind_str = node.annotations.get_value("index")
            node.annotations.add_new(
                name="ancestral",
                value=asrs.loc[i, "end_" + node_ind_str].replace(",", "")
            )

        trees.append(tree)

    # Output the BEAST-formatted trees.
    trees.write_to_path(dest=output_base + "_beast_format.trees", schema="nexus")
