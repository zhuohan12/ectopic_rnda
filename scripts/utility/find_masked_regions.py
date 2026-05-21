#!/usr/bin/env python3
"""
find_masked_regions.py

Scan a FASTA file and output a BED file of contiguous N/n regions
(hard-masked segments).

Usage:
    python find_masked_regions.py input.fa output_masked.bed

Each BED line: chrom  start  end
"""

import sys
import gzip

def open_file(path):
    """Open plain or gzipped FASTA."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def find_masked_regions(fasta_path, bed_path):
    with open_file(fasta_path) as fa, open(bed_path, "w") as out:
        chrom = None
        pos = 0         # 1-based coordinate in FASTA
        in_mask = False
        start = None

        for line in fa:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):  # new chromosome/contig
                if in_mask:  # close out previous masked run
                    out.write(f"{chrom}\t{start-1}\t{pos}\n")
                    in_mask = False
                chrom = line[1:].split()[0]
                pos = 0
                continue

            seq = line.upper()
            for base in seq:
                pos += 1
                if base == "N":
                    if not in_mask:
                        start = pos
                        in_mask = True
                elif in_mask:
                    out.write(f"{chrom}\t{start-1}\t{pos-1}\n")
                    in_mask = False

        # close any trailing mask at EOF
        if in_mask:
            out.write(f"{chrom}\t{start-1}\t{pos}\n")

    print(f"Masked regions written to: {bed_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python find_masked_regions.py input.fa output.bed\n")
    find_masked_regions(sys.argv[1], sys.argv[2])
