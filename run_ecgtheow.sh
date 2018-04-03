#!/bin/bash

set -e

# Parse the command line arguments.
NARGS="$#"
ARGS="$@"

if [ "${NARGS}" -eq 0 ]
  then
    echo "ERROR: No command line arguments were specified."
    exit 1
fi

while [ "${NARGS}" -gt 0 ]
do
  case "$1" in
    --data-dir)
      DATA_DIR="${2%/}"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --sample)
      SAMPLE="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --seed)
      SEED="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --beast-dir)
      BEAST_DIR="${2%/}"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --beagle-dir)
      BEAGLE_DIR="${2%/}"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --nprune)
      NPRUNE="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --mcmc-iter)
      MCMC_ITER="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --mcmc-thin)
      MCMC_THIN="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --mcmc-burnin)
      MCMC_BURNIN="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --asr-nfilters)
      ASR_NFILTERS="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --overwrite)
      OVERWRITE=1
      shift 1
      NARGS=$((NARGS-1))
      ;;
    --run-beast)
      RUN_BEAST=1
      shift 1
      NARGS=$((NARGS-1))
      ;;
    --run-revbayes)
      RUN_REVBAYES=1
      shift 1
      NARGS=$((NARGS-1))
      ;;
    --naive-correction)
      NAIVE_CORRECTION="--naive-correction"
      shift 1
      NARGS=$((NARGS-1))
      ;;
    *)
      echo "ERROR: Please specify valid command line arguments."
      exit 1
      ;;
  esac
done

if [ -z "${NPRUNE}" ]; then NPRUNE=100; fi
if [ -z "${MCMC_ITER}" ]; then MCMC_ITER=10000000; fi
if [ -z "${MCMC_THIN}" ]; then MCMC_THIN=1000; fi
if [ -z "${MCMC_BURNIN}" ]; then MCMC_BURNIN=1000; fi
if [ -z "${ASR_NFILTERS}" ]; then ASR_NFILTERS="50,100"; fi
if [ -z "${OVERWRITE}" ]; then OVERWRITE=0; fi
if [ -z "${RUN_BEAST}" ]; then RUN_BEAST=0; fi
if [ -z "${RUN_REVBAYES}" ]; then RUN_REVBAYES=0; fi
if [ -z "${NAIVE_CORRECTION}" ]; then NAIVE_CORRECTION=""; fi


FAIL=0
if [ -z "${DATA_DIR}" ]; then echo "ERROR: Please specify the '--data-dir' command line argument."; FAIL=1; fi
if [ -z "${SAMPLE}" ]; then echo "ERROR: Please specify the '--sample' command line argument."; FAIL=1; fi
if [ -z "${SEED}" ]; then echo "ERROR: Please specify the '--seed' command line argument."; FAIL=1; fi
if [ -z "${BEAST_DIR}" ]; then echo "ERROR: Please specify the '--beast-dir' command line argument."; FAIL=1; fi
if [ -z "${BEAGLE_DIR}" ]; then echo "ERROR: Please specify the '--beagle-dir' command line argument."; FAIL=1; fi
if [ "${FAIL}" -eq 1 ]; then exit 1; fi


# Check if the output directory already exists.
OUTPUT_DIR=${SEED}_nprune${NPRUNE}_iter${MCMC_ITER}_thin${MCMC_THIN}_burnin${MCMC_BURNIN}
if [ "${OVERWRITE}" -eq 0 -a -d "${OUTPUT_DIR}" ]
then
  echo "ERROR: The output directory already exists! To overwrite the directory, please specify the '--overwrite' command line argument."
  exit 1
fi

# Make the data and runs directories.
mkdir -p ${OUTPUT_DIR}/data ${OUTPUT_DIR}/runs

# Print the command line arguments to a file.
echo ${ARGS} > ${OUTPUT_DIR}/args.log

# Parse the partis YAML info file and get the "healthy" sequences.
export PARTIS=${PWD%/}/lib/cft/partis
python/parse_partis_data.py ${DATA_DIR} --sample ${SAMPLE} --seed ${SEED} --output-dir ${OUTPUT_DIR}

# Generate a tree and prune sequences from the clonal family.
FastTree -nt ${OUTPUT_DIR}/data/healthy_seqs.fasta > ${OUTPUT_DIR}/data/healthy_seqs.tre
lib/cft/bin/prune.py --naive naive --seed ${SEED} ${OUTPUT_DIR}/data/healthy_seqs.tre -n ${NPRUNE} --strategy seed_lineage ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.ids
seqmagick convert --include-from-file ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.ids ${OUTPUT_DIR}/data/healthy_seqs.fasta ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.fasta

# Output the FastTree .PNG tree graphic highlighting the pruned nodes.
python/annotate_fasttree_tree.py ${OUTPUT_DIR}/data/healthy_seqs.tre ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.ids --naive naive --seed ${SEED}

# Trim off site columns with full N-padding.
awk '/^[^>]/ {gsub("N", "-", $0)} {print}' < ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.fasta > ${OUTPUT_DIR}/data/temp.fasta
seqmagick mogrify --squeeze ${OUTPUT_DIR}/data/temp.fasta
awk '/^[^>]/ {gsub("-", "N", $0)} {print}' < ${OUTPUT_DIR}/data/temp.fasta > ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.fasta
rm ${OUTPUT_DIR}/data/temp.fasta

# For both BEAST and RevBayes:
# 1) Construct the input file;
# 2) Run the program;
# 3) Summarize the results.
if [ "${RUN_BEAST}" -eq 1 ]
then
  python/generate_beast_xml_input.py templates/beast_template.xml ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.fasta --naive naive --iter ${MCMC_ITER} --thin ${MCMC_THIN} --output-dir ${OUTPUT_DIR} ${NAIVE_CORRECTION}
  java -Xms64m -Xmx2048m -Djava.library.path=${BEAGLE_DIR} -Dbeast.plugins.dir=beast/plugins -jar ${BEAST_DIR}/lib/beast.jar -warnings -seed 1 -overwrite ${OUTPUT_DIR}/runs/healthy_seqs_nprune${NPRUNE}_beast.xml
  python/trees_to_counted_ancestors.py ${OUTPUT_DIR}/runs/healthy_seqs_nprune${NPRUNE}_beast.trees ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.fasta --naive naive --seed ${SEED} --burnin ${MCMC_BURNIN} --filters ${ASR_NFILTERS//,/ }
fi

if [ "${RUN_REVBAYES}" -eq 1 ]
then
  python/generate_rb_rev_input.py templates/rb_template.rev ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.fasta --naive naive --iter ${MCMC_ITER} --thin ${MCMC_THIN} --output-dir ${OUTPUT_DIR} ${NAIVE_CORRECTION}
  lib/revbayes/projects/cmake/rb ${OUTPUT_DIR}/runs/healthy_seqs_nprune${NPRUNE}_rb.rev
  python/revbayes_to_beast_trees.py ${OUTPUT_DIR}/runs/healthy_seqs_nprune${NPRUNE}_rb.trees ${OUTPUT_DIR}/runs/healthy_seqs_nprune${NPRUNE}_rb.ancestral_states.log
  python/trees_to_counted_ancestors.py ${OUTPUT_DIR}/runs/healthy_seqs_nprune${NPRUNE}_rb_beast.trees ${OUTPUT_DIR}/data/healthy_seqs_nprune${NPRUNE}.fasta --naive naive --seed ${SEED} --burnin ${MCMC_BURNIN} --filters ${ASR_NFILTERS//,/ }
fi
