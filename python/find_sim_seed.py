#!/usr/bin/env python

import argparse
import pickle


def load_tree(filename):
    with open(filename, 'rb') as fh:
        return pickle.load(fh)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('simulated_tree_pickle')
    parser.add_argument('seedid_file')
    return parser.parse_args()


def main():
    args = get_args()
    tree = load_tree(args.simulated_tree_pickle)
    seed, _ = tree.get_farthest_leaf()
    with open(args.seedid_file, 'w') as fh:
        fh.write(seed.name)

if __name__ == '__main__':
    main()


