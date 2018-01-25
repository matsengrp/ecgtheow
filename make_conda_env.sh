#!/bin/bash

# Create the conda environment.
conda create -y -n ecgtheow
source activate ecgtheow
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
printf '#!/bin/sh\n\nexport PYTHONNOUSERSITE=1' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
printf '#!/bin/sh\n\nunset PYTHONNOUSERSITE' > $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

conda install -y biopython graphviz jinja2 psutil python-graphviz pyyaml
conda install -y -c bioconda colorbrewer dendropy seqmagick
conda install -y -c etetoolkit ete3
