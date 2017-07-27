#!/bin/bash

set -e

# Parse the command line arguments.
NARGS="$#"

if [ "${NARGS}" -eq 0 ]
  then
    echo "ERROR: No command line arguments were specified."
    exit 1
fi

while [ "${NARGS}" -gt 0 ]
do
  case "$1" in
    -n|--naive)
      NAIVE="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    -s|--seed)
      SEED="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --beast-dir)
      BEAST_DIR="${2%/}"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --clone-dir)
      CLONE_DIR="${2%/}"
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
    --asr-count-filter)
      ASR_COUNT_FILTER="$2"
      shift 2
      NARGS=$((NARGS-2))
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
if [ -z "${ASR_COUNT_FILTER}" ]; then ASR_COUNT_FILTER=100; fi

if [ -z "${NAIVE}" ] || [ -z "${SEED}" ] || [ -z "${BEAST_DIR}" ] || [ -z "${CLONE_DIR}" ]
  then
    echo "ERROR: Please specify the '--naive', '--seed', '--beast-dir', and '--clone-dir' arguments."
    exit 1
fi

# Make the data and runs directories.
mkdir -p data runs

# Grab the sequences from stoat.
scp stoat:/fh/fast/matsen_e/processed-data/partis/laura-mb/v9/seeds/${SEED}/BF520-h-IgG/run-viterbi-best-plus-0.csv data/${SEED}.csv

# Create the unpruned FASTA file (containing both the naive and seed sequences).
${CLONE_DIR}/pandis/transpose_family.py data/${SEED}.csv
${CLONE_DIR}/pandis/healthy_to_fasta.py data/${SEED}.family_0.csv
${CLONE_DIR}/pandis/get_naives.py data/${SEED}.csv >> data/${SEED}.family_0.healthy.fasta

# Generate a tree and prune sequences from the clonal family.
FastTree -nt data/${SEED}.family_0.healthy.fasta > data/${SEED}.family_0.healthy.tre
${CLONE_DIR}/cft/bin/prune.py --naive ${NAIVE} --seed ${SEED} data/${SEED}.family_0.healthy.tre -n ${NPRUNE} --strategy seed_lineage > data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.ids
seqmagick convert --include-from-file data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.ids data/${SEED}.family_0.healthy.fasta data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta

# Construct the BEAST XML input file.
python/generate_beast_xml_input.py --naive ${NAIVE} --seed ${SEED} templates/beast_template.xml data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta --iter ${MCMC_ITER} --thin ${MCMC_THIN}

# Run BEAST.
java -Xms64m -Xmx2048m -Dbeast.plugins.dir=beast/plugins -jar ${BEAST_DIR}/lib/beast.jar -warnings -seed 1 -overwrite runs/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.xml
if [ -e "${SEED}.family_0.healthy.seedpruned.${NPRUNE}.log" ]
  then
    mv ${SEED}.family_0.healthy.seedpruned.${NPRUNE}.* runs/
fi

# Summarize the BEAST results.
python/trees_to_counted_ancestors.py --seed ${SEED} --burnin ${MCMC_BURNIN} --filter ${ASR_COUNT_FILTER} runs/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.trees data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta
