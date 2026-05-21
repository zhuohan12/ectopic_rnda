import os
import pandas as pd
from collections import defaultdict

def build_blacklist_from_cohort(input_dir, output_prefix="blacklist", 
                                 window_size=500, min_sample_freq=0.1,
                                 file_pattern=".txt"):
    """
    Build a blacklist of recurrent insertion sites across a cohort.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing insertion tables (one per sample)
    output_prefix : str
        Prefix for output files
    window_size : int
        Merge insertions within this distance (bp)
    min_sample_freq : float
        Minimum fraction of samples for blacklisting (0-1)
    file_pattern : str
        File extension pattern to match
    
    Returns:
    --------
    DataFrame with blacklist regions and their frequencies
    """
    
    all_insertions = []
    sample_names = []
    
    for filename in sorted(os.listdir(input_dir)):
        if filename.endswith(file_pattern):
            filepath = os.path.join(input_dir, filename)
            sample_name = filename.replace(file_pattern, "")
            sample_names.append(sample_name)
            
            try:
                df = pd.read_csv(filepath, sep='\t', header=0,
                    usecols=['chrom', 'start', 'end'],
                    dtype={'chrom': str, 'start': int, 'end': int})
                df['sample'] = sample_name
                all_insertions.append(df)
                print(f"Loaded {len(df)} insertions from {filename}")
            except Exception as e:
                print(f"Error reading {filename}: {e}")
    
    if not all_insertions:
        print("No insertion files found!")
        return None
    
    combined = pd.concat(all_insertions, ignore_index=True)
    total_samples = len(sample_names)
    print(f"\nTotal: {len(combined)} insertions from {total_samples} samples")
    
    combined = combined.sort_values(['chrom', 'start', 'end']).reset_index(drop=True)
    
    blacklist_regions = []
    
    for chrom in combined['chrom'].unique():
        chrom_data = combined[combined['chrom'] == chrom].copy()
        
        if len(chrom_data) == 0:
            continue
        
        chrom_data = chrom_data.sort_values('start')
        
        clusters = []
        current_cluster = {
            'chrom': chrom,
            'start': chrom_data.iloc[0]['start'],
            'end': chrom_data.iloc[0]['end'],
            'samples': {chrom_data.iloc[0]['sample']}
        }
        
        for _, row in chrom_data.iloc[1:].iterrows():
            if row['start'] <= current_cluster['end'] + window_size:
                current_cluster['end'] = max(current_cluster['end'], row['end'])
                current_cluster['samples'].add(row['sample'])
            else:
                clusters.append(current_cluster)
                current_cluster = {
                    'chrom': chrom,
                    'start': row['start'],
                    'end': row['end'],
                    'samples': {row['sample']}
                }
        
        clusters.append(current_cluster)
        blacklist_regions.extend(clusters)
    
    blacklist_df = pd.DataFrame([
        {
            'chrom': r['chrom'],
            'start': int(r['start']),
            'end': int(r['end']),
            'width': int(r['end'] - r['start']),
            'n_samples': len(r['samples']),
            'sample_freq': len(r['samples']) / total_samples,
            'samples': ','.join(sorted(r['samples']))
        }
        for r in blacklist_regions
    ])
    
    # Sort by frequency descending — preserves full distribution for plotting
    blacklist_df = blacklist_df.sort_values('n_samples', ascending=False).reset_index(drop=True)
    
    # Summary statistics
    print(f"\n{'='*60}")
    print("BLACKLIST SUMMARY")
    print(f"{'='*60}")
    print(f"Total unique regions: {len(blacklist_df)}")
    print(f"Total samples: {total_samples}")
    print(f"\nFrequency distribution:")
    print(blacklist_df['n_samples'].describe())
    
    print(f"\nRegions by sample count:")
    for threshold in [1, 2, 3, 5, 10, 20]:
        count = (blacklist_df['n_samples'] >= threshold).sum()
        pct = count / len(blacklist_df) * 100 if len(blacklist_df) > 0 else 0
        print(f"  >= {threshold} samples: {count} regions ({pct:.1f}%)")
    
    blacklist_filtered = blacklist_df[blacklist_df['sample_freq'] >= min_sample_freq].copy()
    print(f"\nBlacklist candidates (>= {min_sample_freq*100:.0f}% of samples): {len(blacklist_filtered)} regions")
    
    print(f"\nTop 20 most recurrent sites:")
    print(blacklist_df[['chrom', 'start', 'end', 'width', 'n_samples', 'sample_freq']].head(20).to_string())
    
    return blacklist_df, blacklist_filtered, total_samples


def save_blacklist(blacklist_df, blacklist_filtered, output_prefix, total_samples):
    """
    Save blacklist files.
    
    Output files:
    -------------
    *_frequency_table.csv      : ALL regions sorted by frequency — use for plotting
                                  (histogram of n_samples, scatter of freq vs position, etc.)
    *_blacklist.bed            : BED file of blacklisted regions (bedtools intersect)
    *_blacklist_detailed.tsv   : Blacklist with sample names and metadata
    """
    
    # ── Full frequency table (CSV for plotting) ──────────────────────────────
    # Contains every region with n_samples and sample_freq columns.
    # Suitable for: histogram of recurrence, rank-frequency plot,
    #               per-chromosome barplot, scatter of position vs frequency.
    freq_output = f"{output_prefix}_frequency_table.csv"
    blacklist_df.to_csv(freq_output, index=False)
    print(f"\nSaved full frequency table (for plotting): {freq_output}")
    print(f"  Columns: chrom, start, end, width, n_samples, sample_freq, samples")
    print(f"  Rows: {len(blacklist_df)} (all regions, sorted by n_samples descending)")

    # ── Blacklist BED (for bedtools -v intersection) ──────────────────────────
    bed_output = f"{output_prefix}_blacklist.bed"
    blacklist_filtered[['chrom', 'start', 'end', 'n_samples']].to_csv(
        bed_output, sep='\t', index=False, header=False
    )
    print(f"Saved blacklist BED: {bed_output}")
    
    # ── Detailed blacklist TSV ────────────────────────────────────────────────
    detailed_output = f"{output_prefix}_blacklist_detailed.tsv"
    blacklist_filtered.to_csv(detailed_output, sep='\t', index=False)
    print(f"Saved detailed blacklist: {detailed_output}")
    
    return freq_output, bed_output, detailed_output


print("="*60)
print("ECTOPIC rDNA BLACKLIST BUILDER")
print("="*60)
print("""
Usage:
------
1. Place all sample insertion tables in a directory
2. Run:

   blacklist_df, blacklist_filtered, n_samples = build_blacklist_from_cohort(
       input_dir="/path/to/insertion_tables/",
       output_prefix="rdna_blacklist",
       window_size=500,        # Merge regions within 500bp
       min_sample_freq=0.10,   # Blacklist if in >=10% of samples
       file_pattern=".txt"     # File extension
   )
   
   save_blacklist(blacklist_df, blacklist_filtered, "rdna_blacklist", n_samples)

Plotting the frequency table:
------------------------------
   import pandas as pd
   import matplotlib.pyplot as plt

   df = pd.read_csv("rdna_blacklist_frequency_table.csv")

   # Histogram of recurrence across cohort
   df['n_samples'].plot(kind='hist', bins=30)

   # Rank-frequency plot (log scale)
   df['n_samples'].reset_index(drop=True).plot(logy=True)

   # Per-chromosome frequency
   df.groupby('chrom')['n_samples'].sum().sort_values().plot(kind='barh')

To filter your calls:
---------------------
bedtools intersect -v -a your_calls.bed -b rdna_blacklist_blacklist.bed > filtered_calls.bed
""")

blacklist_df, blacklist_filtered, n_samples = build_blacklist_from_cohort(
    input_dir="/d/hd09/zhuohan/ectopic_rdna_proj/data/output_CCLE_v2",
    output_prefix="rdna_blacklist_v2",
    window_size=500,
    min_sample_freq=0.01,
    file_pattern=".tsv"
)

save_blacklist(blacklist_df, blacklist_filtered, "rdna_blacklist", n_samples)