# Requirements

BEAST (http://beast.community/)

FastTree (http://www.microbesonline.org/fasttree/)

Python modules:

    conda install -y -c etetoolkit ete3 ete3_external_apps
    conda install -y biopython graphviz python-graphviz pandas
    conda install -y -c bioconda colorbrewer dendropy=4.2.0
    pip install jinja2
    pip install seqmagick


# Example usage

Clone the ecgtheow repository using the following command:

    git clone --recursive https://github.com/matsengrp/ecgtheow.git

The `run_beast_asr.sh` shell script in the root directory executes the entire ancestral sequence reconstruction pipeline.  An example is given below:

    ./run_beast_asr.sh --naive naive0 --seed BF520.1-igh --beast-dir /usr/local/BEASTv1.8.4/ --beagle-dir /usr/local/lib/ --data-path /fh/fast/matsen_e/processed-data/partis/laura-mb/v9/seeds/BF520.1-igh/BF520-h-IgG/run-viterbi-best-plus-0.csv

The full set of command line arguments is described in the following table:

| Command | Description |
| ---     | ---         |
| `--naive` | The name of the naive sequence. |
| `--seed` | The name of the seed sequence. |
| `--beast-dir` | The absolute path of the BEAST program. |
| `--beagle-dir` | The absolute path of the BEAGLE libraries. |
| `--data-path` | The absolute path of the partis data file. |
| `--nprune` | The number of sequences to keep from the `cft` pruning step (defaults to 100). |
| `--mcmc-iter` | The number of total MCMC iterations run in BEAST (defaults to 10000000). |
| `--mcmc-thin` | The MCMC sampling frequency used in BEAST (defaults to 1000). |
| `--mcmc-burnin` | The number of MCMC samples thrown away due to burn-in (defaults to 1000). |
| `--asr-nfilter` | The threshold used to filter out infrequent ancestral sequence transitions in the MCMC samples (defaults to 100). |

Note that on stoat, `--beast-dir` and `--beagle-dir` should be set to `/home/matsengrp/local/BEASTv1.8.4/` and `/home/matsengrp/local/lib/`, respectively.
