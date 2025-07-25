#!/bin/bash

# Usage check
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <output_dir> <reference_genome.fa> [target.fa]"
    exit 1
fi

outdir="$1"
background="$2"
target="${3:-/d/hd09/zhuohan/ectopic_rdna_proj/data/U13369.1.fasta}"  # Default to rDNA fasta if not provided

mkdir -p "$outdir"

for k in {17..19..2}
do
    prefix="$(basename "$background" .fa)_comp_k${k}"
    kat comp -m $k -H 10000000000 -t 16 \
        -o "$outdir/$prefix" \
        "$background" "$target"
done
