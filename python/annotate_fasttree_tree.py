#!/usr/bin/env python

import argparse
import os
from ete3 import Tree, TextFace, TreeStyle


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Annotate FastTree tree with prune ids.")
    parser.add_argument(
        'tree_path', type=str,
        help="Path to FastTree tree file.")
    parser.add_argument(
        'ids_path', type=str,
        help="Path to prune id file.")

    args = parser.parse_args()

    tree = Tree(args.tree_path, format=1)

    with open(args.ids_path) as f:
        ids = f.readlines()
    ids = [id.rstrip("\n") for id in ids]

    for leaf_node in tree.get_leaves():
        node_face = TextFace(leaf_node.name, fsize=10,
                             fgcolor="red" if leaf_node.name in ids else "black")
        leaf_node.add_face(node_face, column=0, position="aligned")

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_scale = False

    tree.render(args.ids_path + ".tre.png", tree_style=ts)
