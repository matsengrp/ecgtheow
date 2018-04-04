#!/usr/bin/env python

import argparse
import pickle
import collections as coll
import dendropy as dendro
import csv
import ete3
import itertools as it
from Bio.Seq import Seq


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('true_treefile')
    parser.add_argument('beastfile')
    parser.add_argument('--naive')
    parser.add_argument('--seed')
    parser.add_argument('--burnin', type=int)
    parser.add_argument('samplesout')
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
    return thread_last(lineage_seqs,
            [it.dropwhile, lambda x: x == naive],
            # Can't reverse a generator :(
            [list],
            [reversed],
            [it.dropwhile, lambda x: x == seed],
            [list],
            [reversed])

def translation(seq_str):
    return str(Seq(seq_str).translate())

def node_seq(n, seq_type='dna'):
    dna_seq = n.annotations.get_value('ancestral')
    return translation(dna_seq) if seq_type == 'aa' else dna_seq

def posterior_sample_data(i, sample_lineage, true_internal_seqs, true_naive_seq, seq_type='dna'):
    internal_seqs = set(internal_lineage_seqs([node_seq(n, seq_type=seq_type) for n in sample_lineage], naive=true_naive_seq))
    n_internal_seqs = len(internal_seqs)
    n_correct_seqs = len(internal_seqs.intersection(true_internal_seqs))
    n_incorrect_seqs = n_internal_seqs - n_correct_seqs
    return dict(sample_step=i,
                # this one gets ommitted in csv output, but is needed for intermediate steps
                seq_type=seq_type,
                internal_seqs=internal_seqs,
                n_internal_seqs=n_internal_seqs,
                n_correct_seqs=n_correct_seqs,
                n_missed_seqs=len(true_internal_seqs) - n_correct_seqs,
                n_incorrect_seqs=n_incorrect_seqs,
                root_seq=node_seq(sample_lineage[0], seq_type))


def posterior_results(args, seq_type="dna"):

    # Get the sequences
    true_naive_seq = args.true_tree.sequence
    true_seed_node = args.true_tree.get_leaves_by_name(args.seed)[0]
    true_lineage_seqs = [n.sequence for n in true_seed_node.get_ancestors()]
    true_seed_seq = true_seed_node.sequence

    # map apply aa mappings according to `seq_type`
    if seq_type == 'aa':
        true_lineage_seqs = map(translation, true_lineage_seqs)
        true_seed_seq = translation(true_seed_seq)
        true_naive_seq = translation(true_naive_seq)

    # compute the set of internal sequences
    true_internal_seqs = set(internal_lineage_seqs(true_lineage_seqs, naive=true_naive_seq))

    n_true_internal_seqs = len(true_internal_seqs)

    seed_dist = hammond_dist(true_seed_seq, true_naive_seq)

    # Process posterior tree data into sample dicts
    samples_data = [
            posterior_sample_data(i, sample_lineage, true_internal_seqs, true_naive_seq, seq_type=seq_type)
            for i, sample_lineage in enumerate(tree_lineages(args))
            if i >= args.burnin]

    # post process a bit, counting samples total, and counting seq occurances
    n_trees = len(samples_data)
    counter = coll.Counter()
    for sample_data in samples_data:
        counter.update(sample_data['internal_seqs'])

    # aggregate data per internal sequence obvserved in posterior
    seqs_data = [
        dict(sequence=s,
            seq_type=seq_type,
            times_seen=n,
            prob=float(n) / n_trees,
            correct=(s in true_internal_seqs),
            # This is an internal seq from a posterior tree; if there are no true internal seqs, there's no
            # well defined truth distance, so return None
            truth_dist=min(hammond_dist(s, ts) for ts in true_internal_seqs) if true_internal_seqs else None,
            # These two are really unique per simulated tree, not per seq; but here for convenience
            n_true_internal_seqs=n_true_internal_seqs,
            seed_dist=seed_dist)
        for s, n in counter.most_common()]

    return (samples_data, seqs_data)


def main():
    args = get_args()

    samples_handle = file(args.samplesout, "w")
    samples_writer = csv.DictWriter(samples_handle,
            extrasaction='ignore',
            fieldnames=["seq_type", "sample_step", "n_internal_seqs", "n_correct_seqs", "n_incorrect_seqs", "n_missed_seqs", "root_seq"])

    # Set up outfile handle
    seqs_handle = file(args.seqsout, "w")
    seqs_writer = csv.DictWriter(seqs_handle,
            extrasaction='ignore',
            fieldnames=["seq_type", "times_seen", "prob", "correct", "truth_dist", "n_true_internal_seqs", "seed_dist", "sequence"])
    seqs_writer.writeheader()

    # Get the data out

    # first the true tree
    with open(args.true_treefile) as fh:
        args.true_tree = pickle.load(fh)

    # compute the data, for dna and aa
    dna_samples_data, dna_seqs_data = posterior_results(args, seq_type='dna')
    aa_samples_data, aa_seqs_data   = posterior_results(args, seq_type='aa')

    # write the data
    samples_writer.writerows(it.chain(dna_samples_data, aa_samples_data))
    seqs_writer.writerows(it.chain(dna_seqs_data, aa_seqs_data))

    seqs_handle.close()
    samples_handle.close()


if __name__ == '__main__':
    main()


