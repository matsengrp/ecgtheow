# Requirements

Conda (Miniconda) package manager

C/C++ compiler

Java SE Runtime Environment (version 1.8.0.92 or above)

BEAST (version 1.8.4 or above) (http://beast.community/)

FastTree (http://www.microbesonline.org/fasttree/)

An X server or xvfb-run, necessary for ete3

# Installation

Clone the repository using the following command:

    git clone --recursive https://github.com/matsengrp/ecgtheow.git

Then, create the ecgtheow conda environment and compile RevBayes by running `./INSTALL` in the ecgtheow root directory.

# Example usage

Running `scons` executes the entire ancestral sequence reconstruction pipeline.  An example is given below:

    source activate ecgtheow && cd ecgtheow
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
| `--run-beast` | Should we run BEAST inference? |
| `--run-revbayes` | Should we run RevBayes inference? |
| `--mcmc-iter` | How many RevBayes MCMC iterations (or 100x BEAST iterations) should we use (defaults to 100000)? |
| `--mcmc-thin` | What RevBayes MCMC thinning frequency (or 100x BEAST thinning frequency) should we use (defaults to 10)? |
| `--mcmc-burnin` | What amount of MCMC burnin should we use (defaults to 1000)? |
| `--asr-nfilters` | The comma-separated list of thresholds used to filter out infrequent ancestral sequence transitions in the MCMC samples (defaults to 50,100). |

The entire set of command line arguments can be found by running `scons --help`.
