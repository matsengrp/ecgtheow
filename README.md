This repository contains software components used for the manuscript *Rapid development of an infant-derived HIV-1 broadly neutralizing antibody lineage* by Simonich, Doepker, et al.
It is not configured for use on data besides that used in the paper, although it could be adapted for other data with some work.

# Requirements

* Conda (Miniconda) package manager
* C++ compiler
* Java SE Runtime Environment (version 1.8.0.92 or above)
* BEAST (version 1.8.4 or above) (http://beast.community/)
* FastTree (http://www.microbesonline.org/fasttree/)
* An X server or xvfb-run, necessary for ete3

# Installation

Clone the repository using the following command:

    git clone --recursive https://github.com/matsengrp/ecgtheow.git

Then, create the ecgtheow conda environment and compile RevBayes by running `./INSTALL` in the ecgtheow root directory.

# Example usage

Running `scons` executes the entire ancestral sequence reconstruction pipeline.  An example is given below:

    source activate ecgtheow
    scons --simulate-data --T=50,105,150 --n=5,90,5 --lambda=2 --lambda0=0.365 --target-dist=15 --target-count=30 --carry-cap=1000 --skip-update=100 --nsims=1 --run-beast --run-revbayes

Remember that the conda environment must be activated before running `scons`.
The most important command line arguments are described in the following table:

| Command | Description |
| ---     | ---         |
| `--simulate-data` | Should we generate GC-simulated data? |
| `--nsims` | The number of GC simulation runs (defaults to 1). |
| `--cft-data` | Should we use CFT data? |
| `--data-dir` | For CFT, what data directory should we use? |
| `--sample` | For CFT, what data sample should we use? |
| `--seed` | For CFT, what seed should we use? |
| `--nprune` | For CFT, how many sequences should we keep from the clonal family (defaults to 100)? |
| `--run-beast` | Should we run BEAST inference? |
| `--run-revbayes` | Should we run RevBayes inference? |
| `--mcmc-iter` | How many RevBayes MCMC iterations (100x for BEAST) should we use (defaults to 10000)? |
| `--mcmc-thin` | What RevBayes MCMC thinning frequency (100x for BEAST) should we use (defaults to 10)? |
| `--mcmc-burnin` | What amount of MCMC burnin should we use (defaults to 100)? |
| `--asr-nfilters` | The comma-separated list of thresholds used to filter out infrequent ancestral sequence transitions in the MCMC samples (defaults to 50,100). |

The entire set of command line arguments can be found by running `scons --help`.
