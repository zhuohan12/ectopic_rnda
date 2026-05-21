#!/bin/bash

# Exit on error
set -e

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 < reference.fa>"
    exit 1
fi

ref_fa="$1"
prefix=$(basename "$ref_fa" .fa)

# Index FASTA 
samtools faidx "$ref_fa"

# Create list of chromosome names, excluding chromosome R
cut -f1 "${ref_fa}.fai" | grep -v '^chrR$' > "keep_${prefix}.list"

#Extract all sequences except chromosome R
samtools faidx "$ref_fa" $(cat "keep_${prefix}.list") > "${prefix}_no_rdna.fa"

echo "${prefix}_no_rdna.fa created"