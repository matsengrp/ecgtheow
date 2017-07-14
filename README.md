# Guidelines

Please only commit changes to Jupyter/Sage notebooks after clearing output.


# Requirements

- biopython.org
- https://pythonhosted.org/DendroPy/
- http://graphviz.readthedocs.io/en/stable/index.html


# Example usage

First pull the runs, etc.

    ./pull-ignored.sh

Then we process the BEAST runs into a graph and a list of sequences:

    python/trees_to_counted_ancestors.py --seed BF520.1-igh --burnin 1000 --filter 100 runs/2017-07-10/BF520.1-h-IgH.family_0.healthy.tre.seedpruned.100.ids.trees data/BF520.1-h-IgH.family_0.healthy.tre.seedpruned.100.ids.fasta

(this takes about 10 mins on my chromebook).
