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
        return [{'id': i} for i in range(options['nsims'])]

    @w.add_target()
    def input_seqs(outdir, c):
        sim_setting = c['simulation_setting']
        outbase = path.join(outdir, "sim_seqs")
        sim_seqs = env.Command(
            outbase + ".fasta",
            [sim_setting[x] for x in ['mutability_file', 'substitution_file', 'random_seq_file']],
            "xvfb-run -a lib/bcr-phylo-benchmark/bin/simulator.py --verbose" \
                    +  " --mutability ${SOURCES[0]}" \
                    +  " --substitution ${SOURCES[1]}" \
                    +  " --random_seq ${SOURCES[2]}" \
                    +  " --outbase " + outbase \
                    +  " --lambda " + str(sim_setting['lambda']) \
                    +  " --lambda0 " + str(sim_setting['lambda0']) \
                    +  " --T " + str(sim_setting['T']) \
                    + (" --n " + str(sim_setting['n']) if sim_setting['n'] else "") \
                    +  " --target_dist " + str(sim_setting['target_dist']) \
                    +  " --target_count " + str(sim_setting['target_count']) \
                    +  " --carry_cap " + str(sim_setting['carry_cap']) \
                    +  " --skip_update " + str(sim_setting['skip_update']) \
                    + (" --selection" if sim_setting['selection'] else "") \
                    +  " --random_seed " + str(c['simulation']['id']) \
                    +  " > " + outbase + ".log")
        env.Depends(sim_seqs, "lib/bcr-phylo-benchmark/bin/simulator.py")
        return sim_seqs

elif options["cft_data"]:

    @w.add_target()
    def cluster_seqs(outdir, c):
        outbase = path.join(outdir, "cluster_seqs")
        cluster_seqs = env.Command(
            outbase + ".fasta",
            path.join(options["data_dir"], "info.yaml"),
            "python/parse_partis_data.py " + options["data_dir"] \
                    + " --sample " + options["sample"] \
                    + " --seed " + options["seed"] \
                    + " --output-path " + outbase + ".fasta")
        env.Depends(cluster_seqs, "python/parse_partis_data.py")
        return cluster_seqs

    @w.add_target()
    def input_seqs(outdir, c):
        outbase = path.join(outdir, "trimmed_cluster_seqs")
        trimmed_cluster_seqs = env.Command(
            outbase + ".fasta",
            c["cluster_seqs"],
            "awk \'/^[^>]/ {gsub(\"N\", \"-\", $0)} {print}\' < $SOURCE > temp.fasta;" + \
            "seqmagick mogrify --squeeze temp.fasta;" + \
            "awk \'/^[^>]/ {gsub(\"-\", \"N\", $0)} {print}\' < temp.fasta > $TARGET;" + \
            "rm temp.fasta;")
        return trimmed_cluster_seqs



#@w.add_target()
#def _simulation_posterior_seqs(outdir, c):
#    return dict()
#@w.add_target()
#def _simulation_posterior_samples(outdir, c):
#    return dict()
# Need to add this to our tripl nestly wrapper
#w.add_aggregate('simulation_posterior_seqs', dict)

# Recording software versions
# ---------------------------

#software_versions.add_software_versions(w)
#
        
#
#
#@w.add_target()
#def _process_sim(outdir, c):
#    return env.Command(
#        [path.join(outdir, x) for x in ['sampled_seqs.fasta', 'seedid']],
#        c['tree'],
#        "python/process_sim.py $SOURCE $TARGETS")
#
#
#@w.add_target()
#def sampled_seqs(outdir, c):
#    temp_fasta = path.join(outdir, "temp.fasta")
#    return env.Command(
#        path.join(outdir, 'trimmed_sampled_seqs.fasta'),
#        c['_process_sim'][0],
#        "awk \'/^[^>]/ {gsub(\"N\", \"-\", $0)} {print}\' < $SOURCE > " + temp_fasta + ";" + \
#        "seqmagick mogrify --squeeze " + temp_fasta + ";" + \
#        "awk \'/^[^>]/ {gsub(\"-\", \"N\", $0)} {print}\' < " + temp_fasta + " > $TARGET;" + \
#        "rm " + temp_fasta + ";")
#
#@w.add_target()
#def seed(outdir, c):
#    return c['_process_sim'][1]
#
#
#
#@w.add_nest(full_dump=True)
#def method(c):
#    return [{'id': method + ("-naive-corrected" if naive_correction else ''),
#             'tool': method,
#             'naive_correction': naive_correction}
#            for method in ['revbayes', 'beast']
#            for naive_correction in [True, False]]
#
#
##python/generate_beast_xml_input.py --naive naive --seed ${SEED} templates/beast_template.xml data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta --iter ${MCMC_ITER} --thin ${MCMC_THIN}
#@w.add_target()
#def posterior(outdir, c):
#    tool = c['method']['tool']
#    naive_correction = c['method']['naive_correction']
#    # Set mcmc iters based on whether or not its a test run
#    base_opts = " $SOURCES --naive simcell_1" + (" --naive-correction" if naive_correction else '')
#    if options['test_run']:
#        base_opts += " --iter 10000 --thin 10"
#    else:
#        base_opts += " --iter 10000000 --thin 1000"
#    if tool == 'beast':
#        config_file = env.Command(
#            path.join(outdir, "lineage_reconstruction.xml"),
#            ["templates/beast_template.xml", c['sampled_seqs']],
#            "python/generate_beast_xml_input.py --xml-path $TARGET " + base_opts)
#        return env.SRun(
#            path.join(outdir, "lineage_reconstruction.trees"),
#            config_file,
#            "java -Xms64m -Xmx2048m" + \
#            # Need to abstract over this non-sense with env variables or something
#            " -Djava.library.path=/home/matsengrp/local/lib" + \
#            " -Dbeast.plugins.dir=beast/plugins" + \
#            " -jar /home/matsengrp/local/BEASTv1.8.4/lib/beast.jar" + \
#            " -warnings -seed 1 -overwrite $SOURCE")
#    elif tool == 'revbayes':
#        config_file = env.Command(
#            path.join(outdir, 'lineage_reconstruction.rev'),
#            ['templates/rb_template.rev', c['sampled_seqs']],
#            'python/generate_rb_rev_input.py --rev-path $TARGET ' + base_opts)
#        rb_tgt = env.SRun(
#            path.join(outdir, 'lineage_reconstruction.trees'),
#            config_file,
#            'lib/revbayes/projects/cmake/rb $SOURCE')
#        env.Depends(rb_tgt, c['sampled_seqs'])
#        tgt = env.Command(
#            path.join(outdir, 'lineage_reconstruction_beast.trees'),
#            [path.join(outdir, x) for x in ['lineage_reconstruction.trees', 'lineage_reconstruction.ancestral_states.log']],
#            "python/revbayes_to_beast_trees.py $SOURCES --output-path " + path.join(outdir, 'lineage_reconstruction.trees'))
#        return tgt
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
