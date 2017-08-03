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
    --beagle-dir)
      BEAGLE_DIR="${2%/}"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --cft-version)
      CFT_VERSION="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --chain)
      CHAIN="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    --isotype)
      ISOTYPE="$2"
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
    --asr-nfilter)
      ASR_NFILTER="$2"
      shift 2
      NARGS=$((NARGS-2))
      ;;
    *)
      echo "ERROR: Please specify valid command line arguments."
      exit 1
      ;;
  esac
done

if [ -z "${CFT_VERSION}" ]; then CFT_VERSION="v9"; fi
if [ -z "${CHAIN}" ]; then CHAIN="h"; fi
if [ -z "${ISOTYPE}" ]; then ISOTYPE="G"; fi
if [ -z "${NPRUNE}" ]; then NPRUNE=100; fi
if [ -z "${MCMC_ITER}" ]; then MCMC_ITER=10000000; fi
if [ -z "${MCMC_THIN}" ]; then MCMC_THIN=1000; fi
if [ -z "${MCMC_BURNIN}" ]; then MCMC_BURNIN=1000; fi
if [ -z "${ASR_NFILTER}" ]; then ASR_NFILTER=100; fi

FAIL=0
if [ -z "${NAIVE}" ]; then echo "ERROR: Please specify the '--naive' command line argument."; FAIL=1; fi
if [ -z "${SEED}" ]; then echo "ERROR: Please specify the '--seed' command line argument."; FAIL=1; fi
if [ -z "${BEAST_DIR}" ]; then echo "ERROR: Please specify the '--beast-dir' command line argument."; FAIL=1; fi
if [ -z "${BEAGLE_DIR}" ]; then echo "ERROR: Please specify the '--beagle-dir' command line argument."; FAIL=1; fi
if [ "${FAIL}" -eq 1 ]; then exit 1; fi


# Make the data and runs directories.
mkdir -p data runs

# Grab the sequences from stoat.
scp stoat:/fh/fast/matsen_e/processed-data/partis/laura-mb/v9/seeds/${SEED}-ig${CHAIN}/${SEED%%.*}-${CHAIN}-Ig${ISOTYPE}/run-viterbi-best-plus-0.csv data/${SEED}-ig${CHAIN}.csv

# Create the unpruned FASTA file (containing both the naive and seed sequences).
lib/pandis/transpose_family.py data/${SEED}-ig${CHAIN}.csv
lib/pandis/healthy_to_fasta.py data/${SEED}-ig${CHAIN}.family_0.csv
lib/pandis/get_naives.py data/${SEED}-ig${CHAIN}.csv >> data/${SEED}-ig${CHAIN}.family_0.healthy.fasta

# Generate a tree and prune sequences from the clonal family.
FastTree -nt data/${SEED}-ig${CHAIN}.family_0.healthy.fasta > data/${SEED}-ig${CHAIN}.family_0.healthy.tre
lib/cft/bin/prune.py --naive ${NAIVE} --seed ${SEED}-ig${CHAIN} data/${SEED}-ig${CHAIN}.family_0.healthy.tre -n ${NPRUNE} --strategy seed_lineage > data/${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.ids
seqmagick convert --include-from-file data/${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.ids data/${SEED}-ig${CHAIN}.family_0.healthy.fasta data/${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.fasta

# Construct the BEAST XML input file.
python/generate_beast_xml_input.py --naive ${NAIVE} --seed ${SEED}-ig${CHAIN} templates/beast_template.xml data/${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.fasta --iter ${MCMC_ITER} --thin ${MCMC_THIN}

# Run BEAST.
java -Xms64m -Xmx2048m -Djava.library.path=${BEAGLE_DIR} -Dbeast.plugins.dir=beast/plugins -jar ${BEAST_DIR}/lib/beast.jar -warnings -seed 1 -overwrite runs/${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.xml
if [ -e "${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.log" ]
  then
    mv ${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.* runs/
fi

# Summarize the BEAST results.
python/trees_to_counted_ancestors.py --seed ${SEED}-ig${CHAIN} --burnin ${MCMC_BURNIN} --filter ${ASR_NFILTER} runs/${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.trees data/${SEED}-ig${CHAIN}.family_0.healthy.seedpruned.${NPRUNE}.fasta
