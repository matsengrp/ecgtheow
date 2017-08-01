#!/bin/bash

set -eu

# Process the old inference files.
if [ "$#" -eq 1 ]
then
  FILE_PATHS=${1%/}/*.fasta
  for path in ${FILE_PATHS}
  do
    python/dna_to_aa.py "${path}"
  done
else
  echo "ERROR: Please specify a directory of old inference files as argument."
  exit 1
fi
