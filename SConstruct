#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Run ecgtheow simulation/validation

This SConstruct does the following:

* Simulates data using krdav's bcr-phylo-benchmark simulation script
* Takes the sampled sequences and runs ecgtheow on them
* Takes the ecgtheow posterior and compares the internal node probabilies with those from the given simulation

'''


# Basic imports
from __future__ import print_function
import os
import sys
import subprocess
import datetime
import getpass
from os import path

from SCons.Script import Environment


import sconsutils
sconsutils #lint

# Nestly things
# this temporarily switches between a local checkout and whatever is installed
# Uncomment this line for a local checkout
#sys.path.append(path.join(os.getcwd(), 'nestly'))
import nestly
from nestly import scons as nestly_scons

# Tripl data modelling
# Uncomment this line for a local checkout
#sys.path = [path.join(os.getcwd(), 'tripl')] + sys.path
from tripl import nestly as nestly_tripl


# Set up SCons environment
environ = os.environ.copy()
environ['TMPDIR'] =  '/tmp' # for bcr-phylo-benchmark/simulation

#bcr_phylo_benchmark_dir = '../bcr-phylo-benchmark'


env = Environment(ENV=environ)
# Add stuff to PATH
# Java binary on stoat
env.PrependENVPath('PATH', '/app/easybuild/software/Java/1.8.0_92/bin')
#env.PrependENVPath('PATH', path.join(bcr_phylo_benchmark_dir, 'bin'))
#env.PrependENVPath('PATH', 'bin')


## Arguments

import SCons.Script as Script

# Set up command line arguments/options

# simulation arguments

Script.AddOption('--simulate-data',
        dest="simulate_data",
        action="store_true",
        default=False,
        help="Should we generate GC simulated data?")

Script.AddOption('--lambda',
        dest="lambda",
        type='str',
        default="2.0",
        help='Poisson branching parameter - simulation')

Script.AddOption('--lambda0',
        dest="lambda0",
        type='str',
        default="0.365",
        help='baseline mutation rate - simulation')

Script.AddOption('--T',
        dest="T",
        type='str',
        default="15,30,45",
        help='GC sampling times - simulation')

Script.AddOption('--n',
        dest="n",
        type='int',
        default=None,
        help='How many GC cells are downsampled? - simulation')

Script.AddOption('--target-dist',
        dest="target_dist",
        type='str',
        default="10",
        help='Distance from naive to selection target - simulation')

Script.AddOption('--target-count',
        dest="target_count",
        type='str',
        default="10",
        help='Number of selection targets - simulation')

Script.AddOption('--carry-cap',
        dest="carry_cap",
        type='str',
        default="1000",
        help='The GC carrying capacity - simulation')

Script.AddOption('--skip-update',
        dest="skip_update",
        type='str',
        default="100",
        help='The number of GC cells to skip before update - simulation')

Script.AddOption('--nsims',
        dest="nsims",
        type="int",
        default=1,
        help="The number of GC simulation runs - simulation")

# CFT arguments

Script.AddOption("--cft-data",
        dest="cft_data",
        action="store_true",
        default=False,
        help="Should we use CFT data?")

Script.AddOption("--data-dir",
        dest="data_dir",
        type="str",
        default=None,
        help="What data directory should we use? - CFT")

Script.AddOption("--sample",
        dest="sample",
        type="str",
        default=None,
        help="What data sample should we use? - CFT")

Script.AddOption("--seed",
        dest="seed",
        type="str",
        default=None,
        help="What seed should we use? - CFT")

Script.AddOption('--nprune',
        dest="nprune",
        type='str',
        default="100",
        help='How many sequences should we keep from the clonal family? - CFT')

# BEAST/RevBayes inference arguments

Script.AddOption("--run-beast",
        dest="run_beast",
        action="store_true",
        default=False,
        help="Should we run BEAST inference?")

Script.AddOption("--run-revbayes",
        dest="run_revbayes",
        action="store_true",
        default=False,
        help="Should we run RevBayes inference?")

Script.AddOption("--naive-correction",
        dest="naive_correction",
        action="store_true",
        default=False,
        help="Should we run naive corrected (in addition to regular) Bayesian inference?")

Script.AddOption('--mcmc-iter',
        dest="mcmc_iter",
        type='str',
        default="100000",
        help='How many RevBayes MCMC iterations (or 100x BEAST iterations) should we use?')

Script.AddOption('--mcmc-thin',
        dest="mcmc_thin",
        type='str',
        default="10",
        help='What RevBayes MCMC thinning frequency (or 100x BEAST thinning frequency) should we use?')





Script.AddOption('--test',
        dest='test_run',
        action='store_true',
        default=False,
        help="Setting this flag does a quick test run (low mcmc iter count and sample size)")

Script.AddOption('--outdir',
        dest='outdir',
        metavar='DIR',
        default='output',
        help="Directory in which to output results; defaults to `output`")

Script.AddOption('--lazy-metadata',
        dest='lazy_metadata',
        action='store_true',
        default=False,
        help="""Turns off the default `AlwaysBuild` setting on the various json metadata targets. Without `AlwaysBuild`, metadata
        files may not update properly if any of the SConstruct code changed. However, this can be useful for improving build
        time when debugging, or when only code in scripts has changed. If `--lazy-metadata` is used, it's best to run again
        before using any of the output metadata.json files.""")


def get_options(env):
    # prefer realpath so that running latest vs explicit vN doesn't require rerun; also need for defaults below
    return dict(
        # simulation arguments
        simulate_data = env.GetOption("simulate_data"),
        lambda_ = [float(x) for x in env.GetOption("lambda").split(",")],
        lambda0 = [float(x) for x in env.GetOption("lambda0").split(",")],
        T = [int(x) for x in env.GetOption("T").split(",")],
        n = env.GetOption("n"),
        target_dist = [int(x) for x in env.GetOption("target_dist").split(",")],
        target_count = [int(x) for x in env.GetOption("target_count").split(",")],
        carry_cap = [int(x) for x in env.GetOption("carry_cap").split(",")],
        skip_update = [int(x) for x in env.GetOption("skip_update").split(",")],
        nsims = env.GetOption("nsims"),

        # CFT arguments
        cft_data = env.GetOption("cft_data"),
        data_dir = env.GetOption("data_dir"),
        sample = env.GetOption("sample"),
        seed = env.GetOption("seed"),
        nprune = [int(x) for x in env.GetOption("nprune").split(",")],

        # BEAST/RevBayes inference arguments
        run_beast = env.GetOption("run_beast"),
        run_revbayes = env.GetOption("run_revbayes"),
        naive_correction = env.GetOption("naive_correction"),
        mcmc_iter = [int(x) for x in env.GetOption("mcmc_iter").split(",")],
        mcmc_thin = [int(x) for x in env.GetOption("mcmc_thin").split(",")],

        test_run = env.GetOption('test_run'),
        always_build_metadata = not env.GetOption('lazy_metadata'),
        outdir_base = env.GetOption('outdir'))




# Initialize nestly!
# ==================

# This lets us create parameter nests and run the given pipeline for each combination of nesting parameters.
# It also lets us "pop" out of nesting levels in order to aggregate on nested results.

# Here we initialize nestly, and create an scons wrapper for it

build_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def git(*args):
    return subprocess.check_output(['git'] + list(args))

options = get_options(env)

assert env.GetOption("help") or options["simulate_data"] != options["cft_data"], \
    "Exactly one of '--simulate-data' and '--cft-data' must be specified"
assert env.GetOption("help") or options["simulate_data"] or all(options[x] for x in ["data_dir", "sample", "seed"]), \
    "Please specify the '--data-dir', '--sample', and '--seed' arguments"

nest = nestly.Nest()
w = nestly_scons.SConsWrap(nest, options['outdir_base'], alias_environment=env)
w = nestly_tripl.NestWrap(w,
        name='build',
        # Need to base hashing off of this for optimal incrementalization
        metadata={'id': ('test-' if options['test_run'] else '') + 'ecgtheow-build-' + build_time.replace(' ', '-'),
                  'time': build_time,
                  'command': " ".join(sys.argv),
                  'workdir': os.getcwd(),
                  'user': getpass.getuser(),
                  'commit': git('rev-parse', 'HEAD'),
                  # putting this here these could eventually be nests, 
                  # This will be really cool :-)
                  'diff': git('diff'),
                  'status': git('status', '--porcelain')},
        always_build_metadata=options['always_build_metadata'],
        base_namespace='ecgtheow',
        id_attrs=['ecgtheow.dataset:id', 'ecgtheow.build:id'])




###### STEP 1: Prepare data

if options["simulate_data"]:

    @w.add_nest(full_dump=True)
    def simulation_setting(c):
        # Setting this information as a nest like this leaves us with a data trail in our output metadata, and
        # let's us extend this unary collection later in order to run multiple replicas.
        return [{
            'id': 'lambda' + str(lambda_) + "_lambda0-" + str(lambda0) + \
                  "_targetdist" + str(target_dist) + "_targetcount" + str(target_count) + \
                  "_carrycap" + str(carry_cap) + "_skipupdate" + str(skip_update),
            'model': 'S5F',
            'mutability_file': "lib/bcr-phylo-benchmark/motifs/Mutability_S5F.csv",
            'substitution_file': "lib/bcr-phylo-benchmark/motifs/Substitution_S5F.csv",
            'random_seq_file': "lib/bcr-phylo-benchmark/sequence_data/AbPair_naive_seqs.fa",
            'lambda': lambda_,
            'lambda0': lambda0,
            'T': " ".join(str(T) for T in options["T"]),
            'n': options["n"],
            'selection': True,
            'target_dist': target_dist,
            'target_count': target_count,
            'carry_cap': carry_cap,
            'skip_update': skip_update
                } for lambda_ in options['lambda_']
                  for lambda0 in options['lambda0']
                  for target_dist in options['target_dist']
                  for target_count in options['target_count']
                  for carry_cap in options["carry_cap"]
                  for skip_update in options["skip_update"]]

    @w.add_nest()
    def simulation(c):
        return [{'id': "sim" + str(i), "number": i} for i in range(options['nsims'])]

    @w.add_target()
    def simulation_output(outdir, c):
        sim_setting = c['simulation_setting']
        outbase = path.join(outdir, "sim_seqs")
        sim_outp = env.Command(
            [outbase + x for x in ["_lineage_tree.p", ".fasta"]],
            [sim_setting[x] for x in ['mutability_file', 'substitution_file', 'random_seq_file']],
            "xvfb-run -a lib/bcr-phylo-benchmark/bin/simulator.py --verbose" \
                    +  " --mutability ${SOURCES[0]}" \
                    +  " --substitution ${SOURCES[1]}" \
                    +  " --random_seq ${SOURCES[2]}" \
                    +  " --outbase " + outbase \
                    +  " --lambda " + str(sim_setting['lambda']) \
                    +  " --lambda0 " + str(sim_setting['lambda0']) \
                    +  " --T " + str(sim_setting['T']) \
                    + (" --n " + str(sim_setting['n']) if sim_setting['n'] is not None else "") \
                    +  " --target_dist " + str(sim_setting['target_dist']) \
                    +  " --target_count " + str(sim_setting['target_count']) \
                    +  " --carry_cap " + str(sim_setting['carry_cap']) \
                    +  " --skip_update " + str(sim_setting['skip_update']) \
                    + (" --selection" if sim_setting['selection'] else "") \
                    +  " --random_seed " + str(c['simulation']['number']) \
                    +  " > " + outbase + ".log")
        env.Depends(sim_outp, "lib/bcr-phylo-benchmark/bin/simulator.py")
        return sim_outp

    @w.add_target()
    def simulation_tree(outdir, c):
        return c["simulation_output"][0]

    @w.add_target()
    def input_seqs(outdir, c):
        return c["simulation_output"][1]

    @w.add_target()
    def seed(outdir, c):
        seed = env.Command(
            path.join(outdir, "seed.txt"),
            c["simulation_tree"],
            "python/find_sim_seed.py $SOURCE $TARGET")
        env.Depends(seed, "python/find_sim_seed.py")
        return seed

elif options["cft_data"]:

    @w.add_target()
    def cluster_seqs(outdir, c):
        cluster_seqs = env.Command(
            path.join(outdir, "cluster_seqs.fasta"),
            path.join(options["data_dir"], "info.yaml"),
            "python/parse_partis_data.py " + options["data_dir"] \
                    + " --sample " + options["sample"] \
                    + " --seed " + options["seed"] \
                    + " --output-path $TARGET")
        env.Depends(cluster_seqs, "python/parse_partis_data.py")
        env.Depends(cluster_seqs, "lib/cft/bin/process_partis.py")
        return cluster_seqs

    @w.add_target()
    def cluster_fasttree(outdir, c):
        cluster_fasttree = env.Command(
            path.join(outdir, "cluster_fasttree.tree"),
            c["cluster_seqs"],
            "/home/matsengrp/local/bin/FastTree -nt $SOURCE > $TARGET")
        env.Depends(cluster_fasttree, "/home/matsengrp/local/bin/FastTree")
        return cluster_fasttree

    @w.add_target()
    def seed(outdir, c):
        return env.Command(
            path.join(outdir, "seed.txt"),
            None,
            "echo " + options["seed"] + " > $TARGET")

    @w.add_nest()
    def nprune(c):
        return [{"id": "nprune" + str(nprune),
                 "value": nprune}
                for nprune in options["nprune"]]

    @w.add_target()
    def pruned_cluster_seqids(outdir, c):
        pruned_cluster_seqids = env.Command(
            path.join(outdir, "pruned_cluster_seqids.txt"),
            c["cluster_fasttree"],
            "lib/cft/bin/prune.py --naive naive --seed " + options["seed"] + " $SOURCE -n " + str(c["nprune"]["value"]) + " --strategy seed_lineage $TARGET")
        env.Depends(pruned_cluster_seqids, "lib/cft/bin/prune.py")
        return pruned_cluster_seqids

    @w.add_target()
    def pruned_cluster_seqs(outdir, c):
        return env.Command(
            path.join(outdir, "pruned_cluster_seqs.fasta"),
            [c["pruned_cluster_seqids"], c["cluster_seqs"]],
            "seqmagick convert --include-from-file $SOURCES $TARGET")

    @w.add_target()
    def input_seqs(outdir, c):
        return env.Command(
            path.join(outdir, "input_seqs.fasta"),
            c["pruned_cluster_seqs"],
            "awk \'/^[^>]/ {gsub(\"N\", \"-\", $0)} {print}\' < $SOURCE > temp.fasta;" + \
            "seqmagick mogrify --squeeze temp.fasta;" + \
            "awk \'/^[^>]/ {gsub(\"-\", \"N\", $0)} {print}\' < temp.fasta > $TARGET;" + \
            "rm temp.fasta")

    @w.add_target()
    def pruned_cluster_fasttree_png(outdir, c):
        pruned_cluster_fasttree_png = env.Command(
            path.join(outdir, "pruned_cluster_fasttree.png"),
            [c["cluster_fasttree"], c["pruned_cluster_seqids"]],
            "python/annotate_fasttree_tree.py $SOURCES --naive naive --seed " + options["seed"] + " --output-path $TARGET")
        env.Depends(pruned_cluster_fasttree_png, "python/annotate_fasttree_tree.py")
        return pruned_cluster_fasttree_png





###### STEP 2: Perform BEAST/RevBayes inference

@w.add_nest(full_dump=True)
def inference_setting(c):
    return [{'id': program_name + ("-naive-corrected" if naive_correction else "") + \
                   "_iter" + str(mcmc_iter * (100 if program_name == "beast" else 1)) + \
                   "_thin" + str(mcmc_thin * (100 if program_name == "beast" else 1)),
             'program_name': program_name,
             'naive_correction': naive_correction,
             'mcmc_iter': mcmc_iter * (100 if program_name == "beast" else 1),
             'mcmc_thin': mcmc_thin * (100 if program_name == "beast" else 1)}
            for program_name in ["beast", "revbayes"] if options["run_" + program_name]
            for naive_correction in {options["naive_correction"], False}
            for mcmc_iter in options["mcmc_iter"]
            for mcmc_thin in options["mcmc_thin"]]

@w.add_target()
def naive(outdir, c):
    inf_setting = c["inference_setting"]
    return env.Command(
        path.join(outdir, "naive.txt"),
        None,
        "echo " + \
        ("_naive_" if inf_setting["program_name"] == "revbayes" and
                      not inf_setting["naive_correction"] else "naive") + \
        " > $TARGET")

@w.add_target()
def templater_output(outdir, c):
    inf_setting = c["inference_setting"]

    if inf_setting["program_name"] == "beast":
        templater = "python/generate_beast_xml_input.py"
        template = "templates/beast_template.xml"
        outpath = [path.join(outdir, "beast_run.xml")]
    elif inf_setting["program_name"] == "revbayes":
        templater = "python/generate_rb_rev_input.py"
        template = "templates/rb_template.rev"
        outpath = [path.join(outdir, x) for x in ["revbayes_run.rev", "input_seqs_rb.fasta"]]

    return env.Command(
        outpath,
        [templater, template, c["input_seqs"]],
        "$SOURCES --naive naive" + \
        " --iter " + str(inf_setting["mcmc_iter"]) + \
        " --thin " + str(inf_setting["mcmc_thin"]) + \
       (" --naive-correction" if inf_setting["naive_correction"] else "") + \
        " --output-path ${TARGETS[0]}" + \
       (" --output-fasta ${TARGETS[1]}" if inf_setting["program_name"] == "revbayes" else ""))

@w.add_target()
def input_script(outdir, c):
    return c["templater_output"][0]

@w.add_target()
def inference_output(outdir, c):
    inf_setting = c["inference_setting"]

    if inf_setting["program_name"] == "beast":

        stdout_path = path.join(outdir, "beast_run.stdout.log")
        inf_target = env.Command(
            path.join(outdir, "beast_run.trees"),
            c["input_script"],
            "java -Xms64m -Xmx2048m" + \
            " -Djava.library.path=/home/matsengrp/local/lib/" + \
            " -Dbeast.plugins.dir=beast/plugins" + \
            " -jar /home/matsengrp/local/BEASTv1.8.4/lib/beast.jar" + \
            " -warnings -seed 1 -overwrite $SOURCE" + \
            " > " + stdout_path)
        env.Depends(inf_target, "/home/matsengrp/local/BEASTv1.8.4/lib/beast.jar")
        return inf_target

    elif inf_setting["program_name"] == "revbayes":

        stdout_path = path.join(outdir, "revbayes_run.stdout.log")
        inf_target = env.Command(
            [path.join(outdir, "revbayes_run" + x) for x in [".trees", ".ancestral_states.log"]],
            c["input_script"],
            "lib/revbayes/projects/cmake/rb $SOURCE" + \
            " > " + stdout_path)
        env.Depends(inf_target, "lib/revbayes/projects/cmake/rb")
        postinf_target = env.Command(
            path.join(outdir, "revbayes_run.beast.trees"),
            inf_target,
            "python/revbayes_to_beast_trees.py $SOURCES --output-path $TARGET")
        env.Depends(postinf_target, "python/revbayes_to_beast_trees.py")
        return postinf_target






#@w.add_target()
#def _simulation_posterior_seqs(outdir, c):
#    return dict()
#@w.add_target()
#def _simulation_posterior_samples(outdir, c):
#    return dict()
# Need to add this to our tripl nestly wrapper
#w.add_aggregate('simulation_posterior_seqs', dict)
#
#
#@w.add_target()
#def ecgtheow_counted_ancestors(outdir, c):
#    return env.Command(
#        path.join(outdir, 'ecgtheow_counted_ancestors.dnamap'),
#        [c['seed'], c['posterior'], c['sampled_seqs']],
#        'python/trees_to_counted_ancestors.py ${SOURCES[1]} ${SOURCES[2]} ' \
#            + '--seed `cat $SOURCE` --naive simcell_1 --burnin 100 --filters 50')
#
#@w.add_target()
#def _process_posterior(outdir, c):
#    tgt = env.Command(
#        [path.join(outdir, 'lineage_posterior.' + x) for x in ['posterior_samples.csv', 'posterior_seqs.csv']],
#        [c['seed'], c['tree'], c['posterior'], c['sampled_seqs']],
#        'python/process_beast.py ${SOURCES[1]} ${SOURCES[2]} ' \
#            + '--seed `cat $SOURCE` --naive simcell_1 --burnin 100 $TARGETS')
#    env.Depends(tgt, 'python/process_beast.py')
#    return tgt
#
#@w.add_target()
#def posterior_seqs(outdir, c):
#    tgt = c['_process_posterior'][1]
#    c['_simulation_posterior_seqs'][c['simulation']['id']] = tgt
#    return tgt
#
#@w.add_target()
#def posterior_samples(outdir, c):
#    tgt = c['_process_posterior'][0]
#    c['_simulation_posterior_samples'][c['simulation']['id']] = tgt
#    return tgt
#
#
#w.pop()
#w.pop()
#
#
#@w.add_target()
#def combined_posterior_seqs(outdir, c):
#    seq_files = c['_simulation_posterior_seqs']
#    return env.Command(
#        path.join(outdir, 'combined_posterior_seqs.csv'),
#        seq_files.values(),
#        'csvstack -n simulation -g {} $SOURCES > $TARGET'.format(','.join(str(x) for x in seq_files.keys())))
#
#@w.add_target()
#def combined_posterior_samples(outdir, c):
#    tree_files = c['_simulation_posterior_samples']
#    return env.Command(
#        path.join(outdir, 'combined_posterior_samples.csv'),
#        tree_files.values(),
#        'csvstack -n simulation -g {} $SOURCES > $TARGET'.format(','.join(str(x) for x in tree_files.keys())))
#
#
## trigger final metadata build
#w.pop()
#
#
#
