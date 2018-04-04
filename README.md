# Requirements

Conda (Miniconda) package manager

Java SE Runtime Environment (version 1.8.0.92 or above)

BEAST (http://beast.community/)

FastTree (http://www.microbesonline.org/fasttree/)

An X server or xvfb-run, necessary for ete3

# Installation

Clone the repository using the following command:

    git clone --recursive https://github.com/matsengrp/ecgtheow.git

Then, create the ecgtheow conda environment and compile RevBayes by running `./INSTALL` in the ecgtheow root directory.

# Example usage

The `run_ecgtheow.sh` shell script in the root directory executes the entire ancestral sequence reconstruction pipeline.  An example is given below:

    source activate ecgtheow && cd ecgtheow
    ./run_ecgtheow.sh --data-dir /fh/fast/matsen_e/processed-data/partis/qa255-synth/v17 --sample QA255-g-merged --seed QA255.105-Vh --beast-dir /home/matsengrp/local/BEASTv1.8.4/ --beagle-dir /home/matsengrp/local/lib/ --run-beast

Remember that the conda environment must be activated before running the `run_ecgtheow.sh` script.
The full set of command line arguments is described in the following table:

| Command | Description |
| ---     | ---         |
| `--data-dir` | The absolute path of the CFT dataset YAML file. |
| `--sample` | The name of the sample. |
| `--seed` | The name of the seed sequence. |
| `--beast-dir` | The absolute path of the BEAST program. |
| `--beagle-dir` | The absolute path of the BEAGLE libraries. |
| `--nprune` | The number of sequences to keep from the CFT pruning step (defaults to 100). |
| `--mcmc-iter` | The number of total MCMC iterations run in BEAST (defaults to 10000000). |
| `--mcmc-thin` | The MCMC sampling frequency used in BEAST (defaults to 1000). |
| `--mcmc-burnin` | The number of MCMC samples thrown away due to burn-in (defaults to 1000). |
| `--asr-nfilters` | The comma-separated list of thresholds used to filter out infrequent ancestral sequence transitions in the MCMC samples (defaults to 50,100). |
| `--overwrite` | A binary flag that indicates whether to overwrite already existing results. |
| `--run-beast` | A binary flag that indicates whether to run BEAST. |
| `--run-revbayes` | A binary flag that indicates whether to run RevBayes. |
| `--naive-correction` | Should we apply the naive sequence likelihood correction? |

## X Servers & ETE3

Note that ETE3 (as mentioned above) requires an active X-server connection.
If the environment you're running from does not have one, you can simulate one by running with `xvfb-run`.
Doing this for the example above would look like:

    xvfb-run ./run_ecgtheow.sh --data-dir /fh/fast/matsen_e/processed-data/partis/qa255-synth/v17 --sample QA255-g-merged --seed QA255.105-Vh --beast-dir /home/matsengrp/local/BEASTv1.8.4/ --beagle-dir /home/matsengrp/local/lib/
