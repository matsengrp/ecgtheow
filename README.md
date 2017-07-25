# Guidelines

Please only commit changes to Jupyter/Sage notebooks after clearing output.


# Requirements

    conda install -y -c etetoolkit ete3 ete3_external_apps
    conda install -y biopython graphviz python-graphviz pandas
    conda install -y -c bioconda colorbrewer dendropy=4.2.0
    pip install seqmagick


# Example usage

To generate a BEAST XML input file, run the following command:

    python/generate_beast_xml_input.py --naive naive0 --seed BF520.1-igh templates/beast_template.xml data/2017-07-10/BF520.1-h-IgH.family_0.healthy.tre.seedpruned.100.ids.fasta

Now, run BEAST using the previously created XML file:

    cp -r beast/plugins ${BEAST_ROOT}
    mv runs/2017-07-10/BF520.1-h-IgH.family_0.healthy.tre.seedpruned.100.ids.xml ${BEAST_ROOT}
    cd ${BEAST_ROOT}
    java -jar lib/beast.jar -seed 0 BF520.1-h-IgH.family_0.healthy.tre.seedpruned.100.ids.xml
    mv BF520.1-h-IgH.family_0.healthy.tre.seedpruned.100.ids.* ${ECGTHEOW_ROOT}/runs/2017-07-10
    cd ${ECGTHEOW_ROOT}

Alternatively, you can pull the runs from `stoat`:

    ./pull-ignored.sh

Then, we process the BEAST runs into a graph and a list of sequences:

    python/trees_to_counted_ancestors.py --seed BF520.1-igh --burnin 1000 --filter 100 runs/2017-07-10/BF520.1-h-IgH.family_0.healthy.tre.seedpruned.100.ids.trees data/2017-07-10/BF520.1-h-IgH.family_0.healthy.tre.seedpruned.100.ids.fasta

(this takes about 10 mins on my chromebook).
