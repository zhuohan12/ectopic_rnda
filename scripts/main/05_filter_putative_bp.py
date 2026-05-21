#!/usr/bin/env python3
"""
Filter ectopic rDNA insertion calls based on read evidence quality.

Usage:
    python filter_rdna_insertions.py <input.tsv> <output.tsv>

Filter criteria:
    1. chrom != 'chrR'       : Exclude rDNA array itself (not ectopic)
    2. nA + nB >= 3          : At least 3 soft-clip reads (Cases A + B)
    3. nC >= 3               : At least 3 discordant reads (Case C)
                              (ensures independent paired-end support)
    4. mean_mapq >= 20       : Minimum mapping quality threshold
    5. chrR_C_span <= 700    : Range of Case C rDNA positions <= 700bp
                              (removes multi-mapping artefacts where one
                               or more reads map to a distant rDNA locus;
                               span directly detects single-read outliers
                               that MAD would miss by design; threshold
                               derived empirically from long-read
                               validation cohort: 100% specificity)
"""

import sys
import pandas as pd
import numpy as np
from collections import Counter


def filter_insertions(input_file, output_file):
    print(f"Reading {input_file}...")
    df = pd.read_csv(input_file, sep='\t')
    print(f"Loaded {len(df)} candidate insertions")

    # Convert numeric columns (handle '.' as missing)
    non_numeric = ['chrom', 'strand_A', 'strand_B', 'strand_C']
    for col in df.columns:
        if col not in non_numeric:
            df[col] = pd.to_numeric(df[col].replace('.', np.nan), errors='coerce')

    # Derived column
    df['nAB'] = df['nA'] + df['nB']

    # Apply filters
    print("\nApplying filters:")
    print("  (1) chrom != 'chrR'")
    print("  (2) nA + nB >= 3")
    print("  (3) nC >= 3")
    print("  (4) mean_mapq >= 20")
    print("  (5) chrR_C_span <= 700")

    mask_not_chrR = df['chrom'] != 'chrR'
    mask_nAB      = df['nAB'] >= 3
    mask_nC       = df['nC'] >= 3
    mask_mapq     = df['mean_mapq'] >= 20
    mask_span     = df['chrR_C_span'].fillna(np.inf) <= 700

    print(f"\nFilter results (independent, against full input):")
    print(f"  chrom != 'chrR':     {mask_not_chrR.sum():>5}/{len(df)} pass ({mask_not_chrR.sum()/len(df)*100:.1f}%)")
    print(f"  nA + nB >= 3:        {mask_nAB.sum():>5}/{len(df)} pass ({mask_nAB.sum()/len(df)*100:.1f}%)")
    print(f"  nC >= 3:             {mask_nC.sum():>5}/{len(df)} pass ({mask_nC.sum()/len(df)*100:.1f}%)")
    print(f"  mean_mapq >= 20:     {mask_mapq.sum():>5}/{len(df)} pass ({mask_mapq.sum()/len(df)*100:.1f}%)")
    print(f"  chrR_C_span <= 700:  {mask_span.sum():>5}/{len(df)} pass ({mask_span.sum()/len(df)*100:.1f}%)")

    keep = mask_not_chrR & mask_nAB & mask_nC & mask_mapq & mask_span

    print(f"\n{'='*60}")
    print(f"FINAL: {keep.sum()}/{len(df)} insertions pass all filters ({keep.sum()/len(df)*100:.1f}%)")
    print(f"{'='*60}")

    df_filtered = df[keep].copy()
    df_filtered.to_csv(output_file, sep='\t', index=False)
    print(f"\nFiltered insertions written to {output_file}")

    # Summary statistics
    if len(df_filtered) > 0:
        print(f"\nSummary of filtered insertions:")
        print(f"  nAB range:          {df_filtered['nAB'].min():.0f} - {df_filtered['nAB'].max():.0f}")
        print(f"  nC range:           {df_filtered['nC'].min():.0f} - {df_filtered['nC'].max():.0f}")
        print(f"  chrR_C_span range:  {df_filtered['chrR_C_span'].min():.1f} - {df_filtered['chrR_C_span'].max():.1f}")
        print(f"  MAPQ range:         {df_filtered['mean_mapq'].min():.1f} - {df_filtered['mean_mapq'].max():.1f}")
        print(f"  Median MAPQ:        {df_filtered['mean_mapq'].median():.1f}")

        chr_counts = df_filtered['chrom'].value_counts().head(10)
        print(f"\n  Top chromosomes:")
        for chrom, count in chr_counts.items():
            print(f"    {chrom}: {count}")

    # Rejection breakdown - first failing filter per locus
    print(f"\nRejection breakdown (first failing filter per locus):")
    rejected = df[~keep].copy()
    reasons = []
    for _, row in rejected.iterrows():
        if row['chrom'] == 'chrR':
            reasons.append('chrom == chrR')
        elif row['nAB'] < 3:
            reasons.append('nAB < 3')
        elif pd.isna(row['nC']) or row['nC'] < 3:
            reasons.append('nC < 3')
        elif pd.isna(row['mean_mapq']) or row['mean_mapq'] < 20:
            reasons.append('mean_mapq < 20')
        elif pd.isna(row['chrR_C_span']) or row['chrR_C_span'] > 700:
            reasons.append('chrR_C_span > 700')
        else:
            reasons.append('unknown')

    for reason, count in Counter(reasons).most_common():
        print(f"  {reason}: {count} loci")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    filter_insertions(sys.argv[1], sys.argv[2])