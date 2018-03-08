#!/usr/bin/env python

import argparse
import pickle
import collections as coll
import dendropy as dendro
import csv
import ete3
import itertools as it


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('true_treefile')
    parser.add_argument('beastfile')
    parser.add_argument('--naive')
    parser.add_argument('--seed')
    parser.add_argument('--burnin', type=int)
    parser.add_argument('treesout')
    parser.add_argument('seqsout')
    return parser.parse_args()


def node_ancestors(node, descendants=None):
    return (node_ancestors(node.parent_node) if node.parent_node else []) + [node] + (descendants or [])

def tree_lineages(args):
    for tree in dendro.Tree.yield_from_files(files=[args.beastfile], schema='nexus', preserve_underscores=True):
        yield node_ancestors(tree.find_node_with_taxon_label(args.seed))

def hammond_dist(s1, s2):
    return sum(x != y for x, y in zip(s1, s2))


def thread_last(init, *computations):
    """Beautiful, data-driven, functional goodness; Runs the init value through each of the computations in
    turn, as lists of [fn, *args], passing the result through as the last arg to each function call."""
    if computations:
        computation = computations[0]
        fn = computation[0]
        args = list(computation[1:])
        args.append(init)
        return thread_last(fn(*args), *computations[1:])
    else:
        return init

def internal_lineage_seqs(lineage_seqs, naive=None):
    "lineage_seqs should be a list, naive, optional; defaults to first in lineage"
    naive = naive or lineage_seqs[0]
    seed = lineage_seqs[-1]
    #def debug(arg, *args, **kw_args):
        #print "arg:", arg
        #print "args:", args
        #print "kw_args:", kw_args
        #return arg
    return thread_last(lineage_seqs,
            [it.dropwhile, lambda x: x == naive],
            [list],
            [reversed],
            [it.dropwhile, lambda x: x == seed],
            [list],
            [reversed])


def main():
    args = get_args()

    counter = coll.Counter()

    with open(args.true_treefile) as fh:
        true_tree = pickle.load(fh)
    true_naive_seq = true_tree.sequence
    seed_node = true_tree.get_leaves_by_name(args.seed)[0]
    true_lineage_seqs = [n.sequence for n in seed_node.get_ancestors()]

    true_internal_seqs = set(internal_lineage_seqs(true_lineage_seqs, naive=true_naive_seq))
    n_true_internal_seqs = len(true_internal_seqs)

    seed_dist = hammond_dist(seed_node.sequence, true_naive_seq)

    trees_handle = file(args.treesout, "w")
    trees_writer = csv.DictWriter(trees_handle,
            fieldnames=["sample_step", "n_internal_seqs", "n_correct_seqs", "n_incorrect_seqs", "n_missed_seqs", "root_seq"])
    trees_writer.writeheader()

    seq = lambda n: n.annotations.get_value('ancestral')
    for i, tree_lineage in enumerate(tree_lineages(args)):
        if i >= args.burnin:
            internal_seqs = set(internal_lineage_seqs([seq(n) for n in tree_lineage], naive=true_naive_seq))
            n_internal_seqs = len(internal_seqs)
            n_correct_seqs = len(internal_seqs.intersection(true_internal_seqs))
            n_incorrect_seqs = n_internal_seqs - n_correct_seqs
            counter.update(internal_seqs)
            trees_writer.writerow(dict(
                sample_step=i,
                n_internal_seqs=n_internal_seqs,
                n_correct_seqs=n_correct_seqs,
                n_missed_seqs=n_true_internal_seqs - n_correct_seqs,
                n_incorrect_seqs=n_incorrect_seqs,
                root_seq=seq(tree_lineage[0])))

    # Set up outfile handle
    seqs_handle = file(args.seqsout, "w")
    seqs_writer = csv.DictWriter(seqs_handle,
            fieldnames=["sequence", "times_seen", "prob", "correct", "truth_dist", "n_true_internal_seqs", "seed_dist"])
    seqs_writer.writeheader()

    n_trees = i + 1
    seqs_writer.writerows(
        dict(sequence=s,
            times_seen=n,
            prob=float(n) / n_trees,
            correct=(s in true_internal_seqs),
            truth_dist=min(hammond_dist(s, ts) for ts in true_internal_seqs),
            # These two are really unique per simulated tree, not per seq; but here for convenience
            n_true_internal_seqs=n_true_internal_seqs,
            seed_dist=seed_dist)
        for s, n in counter.most_common())

    seqs_handle.close()
    trees_handle.close()


if __name__ == '__main__':
    main()


