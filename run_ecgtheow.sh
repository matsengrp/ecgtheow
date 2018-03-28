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
mkdir -p data runs

# Print the command line call to a file.
echo ${ARGS} > args.log

# Parse the partis YAML info file and get the "healthy" sequences.
export PARTIS=${PWD%/}/lib/cft/partis
python/parse_partis_data.py ${DATA_DIR} --sample ${SAMPLE} --seed ${SEED}

# Generate a tree and prune sequences from the clonal family.
FastTree -nt data/${SEED}.family_0.healthy.fasta > data/${SEED}.family_0.healthy.tre
lib/cft/bin/prune.py --naive naive --seed ${SEED} data/${SEED}.family_0.healthy.tre -n ${NPRUNE} --strategy seed_lineage data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.ids
seqmagick convert --include-from-file data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.ids data/${SEED}.family_0.healthy.fasta data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta

# Output the FastTree .PNG tree graphic highlighting the pruned nodes.
python/annotate_fasttree_tree.py data/${SEED}.family_0.healthy.tre data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.ids --naive naive --seed ${SEED}

# Trim off site columns with full N-padding.
awk '/^[^>]/ {gsub("N", "-", $0)} {print}' < data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta > data/temp.fasta
seqmagick mogrify --squeeze data/temp.fasta
awk '/^[^>]/ {gsub("-", "N", $0)} {print}' < data/temp.fasta > data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta
rm data/temp.fasta

# Construct the BEAST XML input file.
python/generate_beast_xml_input.py --naive naive --seed ${SEED} templates/beast_template.xml data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta --iter ${MCMC_ITER} --thin ${MCMC_THIN}

# Run BEAST.
java -Xms64m -Xmx2048m -Djava.library.path=${BEAGLE_DIR} -Dbeast.plugins.dir=beast/plugins -jar ${BEAST_DIR}/lib/beast.jar -warnings -seed 1 -overwrite runs/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.xml
if [ -e "${SEED}.family_0.healthy.seedpruned.${NPRUNE}.log" ]
then
  mv ${SEED}.family_0.healthy.seedpruned.${NPRUNE}.* runs/
fi

# Summarize the BEAST results.
python/trees_to_counted_ancestors.py runs/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.trees data/${SEED}.family_0.healthy.seedpruned.${NPRUNE}.fasta --seed ${SEED} --burnin ${MCMC_BURNIN} --filters ${ASR_NFILTERS//,/ }

# Move the results to the output directory.
if [ "${OVERWRITE}" -eq 1 ]
then
  mkdir -p ${OUTPUT_DIR}
  cp -t ${OUTPUT_DIR} -r data/ runs/ args.log
  rm -r data/ runs/ args.log
else
  mkdir ${OUTPUT_DIR}
  mv -t ${OUTPUT_DIR} data/ runs/ args.log
fi
