import os
import pandas as pd

def load_blacklist(bed_file):
    """
    Load blacklist regions from BED file.
    
    Parameters:
    -----------
    bed_file : str
        Path to BED file (chrom, start, end, ...)
    
    Returns:
    --------
    dict : {chrom: [(start, end), ...]}
    """
    blacklist = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            if chrom not in blacklist:
                blacklist[chrom] = []
            blacklist[chrom].append((start, end))
    
    # Sort intervals for each chromosome
    for chrom in blacklist:
        blacklist[chrom].sort()
    
    return blacklist


def overlaps_blacklist(chrom, start, end, blacklist):
    """
    Check if a region overlaps any blacklist region.
    
    Parameters:
    -----------
    chrom : str
    start : int
    end : int
    blacklist : dict
    
    Returns:
    --------
    bool : True if overlaps blacklist
    """
    if chrom not in blacklist:
        return False
    
    for bl_start, bl_end in blacklist[chrom]:
        # Check for overlap
        if start <= bl_end and end >= bl_start:
            return True
        # Early exit if past the region (since sorted)
        if bl_start > end:
            break
    
    return False


def filter_insertions(input_dir, output_dir, blacklist_bed, file_pattern=".tsv"):
    """
    Filter insertion TSVs by removing blacklisted regions.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing input TSV files
    output_dir : str
        Directory for filtered output TSV files
    blacklist_bed : str
        Path to blacklist BED file
    file_pattern : str
        File extension to match
    """
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load blacklist
    print(f"Loading blacklist from: {blacklist_bed}")
    blacklist = load_blacklist(blacklist_bed)
    total_bl_regions = sum(len(v) for v in blacklist.values())
    print(f"Loaded {total_bl_regions} blacklist regions across {len(blacklist)} chromosomes")
    
    # Process each file
    summary = []
    
    for filename in sorted(os.listdir(input_dir)):
        if not filename.endswith(file_pattern):
            continue
        
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)
        
        # Read TSV with header
        df = pd.read_csv(input_path, sep='\t', header=0)
        original_count = len(df)
        
        # Filter out blacklisted regions
        mask = df.apply(
            lambda row: not overlaps_blacklist(
                str(row['chrom']), 
                int(row['start']), 
                int(row['end']), 
                blacklist
            ), 
            axis=1
        )
        
        df_filtered = df[mask]
        filtered_count = len(df_filtered)
        removed_count = original_count - filtered_count
        
        # Save filtered TSV (preserving all columns and header)
        df_filtered.to_csv(output_path, sep='\t', index=False)
        
        summary.append({
            'sample': filename,
            'original': original_count,
            'filtered': filtered_count,
            'removed': removed_count,
            'pct_removed': f"{removed_count/original_count*100:.1f}%" if original_count > 0 else "0%"
        })
        
        print(f"{filename}: {original_count} -> {filtered_count} ({removed_count} removed)")
    
    # Print summary
    print(f"\n{'='*60}")
    print("FILTERING SUMMARY")
    print(f"{'='*60}")
    summary_df = pd.DataFrame(summary)
    print(summary_df.to_string(index=False))
    
    total_original = summary_df['original'].sum()
    total_filtered = summary_df['filtered'].sum()
    total_removed = summary_df['removed'].sum()
    print(f"\nTotal: {total_original} -> {total_filtered} ({total_removed} removed, {total_removed/total_original*100:.1f}%)")
    print(f"\nFiltered files saved to: {output_dir}")
    
    return summary_df


# Usage example
print("""
BLACKLIST FILTER FOR ECTOPIC rDNA INSERTIONS
=============================================

Usage:
------
filter_insertions(
    input_dir="/path/to/insertion_tsvs/",
    output_dir="/path/to/filtered_output/",
    blacklist_bed="/path/to/blacklist.bed",
    file_pattern=".tsv"
)

Example:
--------
filter_insertions(
    input_dir="/d/hd09/zhuohan/ectopic_rdna_proj/data/output_CCLE/",
    output_dir="/d/hd09/zhuohan/ectopic_rdna_proj/data/output_CCLE_filtered/",
    blacklist_bed="/d/hd09/zhuohan/ectopic_rdna_proj/data/rdna_blacklist_blacklist.bed",
    file_pattern=".tsv"
)
""")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--blacklist_bed", required=True)
    parser.add_argument("--file_pattern", default=".tsv")
    args = parser.parse_args()

    filter_insertions(args.input_dir, args.output_dir, args.blacklist_bed, args.file_pattern)