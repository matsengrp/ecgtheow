#!/usr/bin/env python

import argparse
import collections
import csv
import dendropy
import ete3
import itertools
import pickle
from util_functions import seqs_of_tree, translate


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('true_treefile')
    parser.add_argument('beastfile')
    parser.add_argument('--seed')
    parser.add_argument('--burnin', type=int)
    parser.add_argument('samplesout')
    parser.add_argument('aggregatesout')
    return parser.parse_args()

def hamming_dist(s1, s2):
    return sum(x != y for x, y in zip(s1, s2))

def process_single_asr(i, asr_lineage_seqs, true_lineage_seqs, seq_type):
    # map apply aa mappings according to `seq_type`
    if seq_type == 'aa':
        asr_lineage_seqs = map(translate, asr_lineage_seqs)

    # keep only the unique lineage sequences
    asr_lineage_seqs = set(asr_lineage_seqs)

    n_correct_seqs = len(asr_lineage_seqs & true_lineage_seqs)
    n_incorrect_seqs = len(asr_lineage_seqs) - n_correct_seqs
    n_missed_seqs = len(true_lineage_seqs) - n_correct_seqs

    return dict(seq_type=seq_type,
                sample_step=i,
                asr_lineage_seqs=asr_lineage_seqs,
                n_correct_seqs=n_correct_seqs,
                n_incorrect_seqs=n_incorrect_seqs,
                n_missed_seqs=n_missed_seqs)

def asr_results(args, seq_type):
    # Get the sequences
    true_naive_seq = args.true_tree.sequence
    true_seed_node = args.true_tree.get_leaves_by_name(args.seed)[0]
    true_lineage_seqs = [n.sequence for n in true_seed_node.iter_ancestors()]
    true_seed_seq = true_seed_node.sequence
    true_lineage_seqs.insert(0, true_seed_seq)

    # map apply aa mappings according to `seq_type`
    if seq_type == 'aa':
        true_lineage_seqs = map(translate, true_lineage_seqs)
        true_seed_seq = translate(true_seed_seq)
        true_naive_seq = translate(true_naive_seq)

    # keep only the unique lineage sequences
    true_lineage_seqs = set(true_lineage_seqs)

    # compute the true naive-seed distance
    true_naive_seed_dist = hamming_dist(true_naive_seq, true_seed_seq)

    # Process posterior tree data into sample dicts
    tree_yielder = dendropy.Tree.yield_from_files(
        files=[args.beastfile],
        schema='nexus',
        preserve_underscores=True
    )
    asr_samples = [
        process_single_asr(i, seqs_of_tree(tree, args.seed) + [args.true_tree.sequence],
                           true_lineage_seqs, seq_type)
        for i, tree in enumerate(tree_yielder) if i >= args.burnin]

    # post process a bit, counting samples total, and counting seq occurances
    n_trees = len(asr_samples)
    counter = collections.Counter()
    for asr_sample in asr_samples:
        counter.update(asr_sample['asr_lineage_seqs'])

    # add missed lineage sequences to the counter
    for seq in true_lineage_seqs:
        if not counter[seq]:
            counter[seq] = 0

    # aggregate data per internal sequence observed in posterior
    asr_aggregates = [
        dict(seq_type=seq_type,
             count=cnt,
             prob=float(cnt) / n_trees,
             correct=(seq in true_lineage_seqs),
             truth_dist=min(hamming_dist(seq, true_seq) for true_seq in true_lineage_seqs),
             # These two are really unique per simulated tree, not per seq; but here for convenience
             n_true_lineage_seqs=len(true_lineage_seqs),
             true_naive_seed_dist=true_naive_seed_dist,
             sequence=seq)
        for seq, cnt in counter.most_common(None)]

    return (asr_samples, asr_aggregates)


def main():
    args = get_args()

    # Set up outfile handles
    samples_handle = file(args.samplesout, "w")
    samples_writer = csv.DictWriter(samples_handle,
                                    extrasaction='ignore',
                                    fieldnames=["seq_type", "sample_step", "n_correct_seqs",
                                                "n_incorrect_seqs", "n_missed_seqs"])
    samples_writer.writeheader()

    aggregates_handle = file(args.aggregatesout, "w")
    aggregates_writer = csv.DictWriter(aggregates_handle,
                                       extrasaction='ignore',
                                       fieldnames=["seq_type", "count", "prob", "correct",
                                                   "truth_dist", "n_true_lineage_seqs",
                                                   "true_naive_seed_dist", "sequence"])
    aggregates_writer.writeheader()

    # Get the data out

    # first the true tree
    with open(args.true_treefile) as fh:
        args.true_tree = pickle.load(fh)

    # compute the data, for dna and aa
    dna_asr_samples, dna_asr_aggregates = asr_results(args, seq_type='dna')
    aa_asr_samples, aa_asr_aggregates   = asr_results(args, seq_type='aa')

    # write the data
    samples_writer.writerows(itertools.chain(dna_asr_samples, aa_asr_samples))
    aggregates_writer.writerows(itertools.chain(dna_asr_aggregates, aa_asr_aggregates))

    # close the file handles
    samples_handle.close()
    aggregates_handle.close()


if __name__ == '__main__':
    main()
