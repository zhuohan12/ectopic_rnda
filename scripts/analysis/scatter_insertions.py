#!/usr/bin/env python3
"""
Enhanced script to analyze ectopic rDNA insertion clusters with multiple visualization options.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
from pathlib import Path
import sys
from scipy.stats import gaussian_kde
import seaborn as sns


def parse_numeric(value):
    """Convert value to numeric, treating '.' as NaN"""
    if value == '.' or value == '' or pd.isna(value):
        return np.nan
    try:
        return float(value)
    except (ValueError, TypeError):
        return np.nan


def process_file(filepath):
    """Process a single TSV file and extract relevant data."""
    try:
        df = pd.read_csv(filepath, sep='\t')
        results = []
        filename = filepath.name
        
        for idx, row in df.iterrows():
            mapq = parse_numeric(row['mean_mapq'])
            if pd.isna(mapq) or mapq <= 20:
                continue
            
            chrR_A_n = parse_numeric(row['chrR_A_n'])
            chrR_B_n = parse_numeric(row['chrR_B_n'])
            chrR_A_span = parse_numeric(row['chrR_A_span'])
            chrR_B_span = parse_numeric(row['chrR_B_span'])
            
            chrR_A_n = 0 if pd.isna(chrR_A_n) else chrR_A_n
            chrR_B_n = 0 if pd.isna(chrR_B_n) else chrR_B_n
            
            total_reads = chrR_A_n + chrR_B_n
            if total_reads < 2:
                continue
            
            spans = []
            if not pd.isna(chrR_A_span):
                spans.append(chrR_A_span)
            if not pd.isna(chrR_B_span):
                spans.append(chrR_B_span)
            
            if len(spans) > 0:
                max_span = max(spans)
                results.append((max_span, total_reads, filename))
        
        return results
    
    except Exception as e:
        print(f"Error processing {filepath}: {e}", file=sys.stderr)
        return []


def create_marginal_plot(total_reads, max_spans, output_path):
    """Create scatterplot with marginal distributions and log-scaled x-axis."""
    
    fig = plt.figure(figsize=(12, 10))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 4, 0.2], width_ratios=[4, 1, 0.2],
                          hspace=0.05, wspace=0.05)
    
    # Main scatterplot
    ax_main = fig.add_subplot(gs[1, 0])
    
    # Calculate point density for coloring
    # Use log scale for both x and y
    log_y = np.log10(np.maximum(max_spans, 0.1))  # Handle zeros
    xy = np.vstack([np.log10(total_reads), log_y])
    z = gaussian_kde(xy)(xy)
    
    idx = z.argsort()
    x_sorted = np.array(total_reads)[idx]
    y_sorted = np.array(max_spans)[idx]
    z_sorted = z[idx]
    
    scatter = ax_main.scatter(x_sorted, y_sorted, c=z_sorted, s=50, 
                             cmap='viridis', alpha=0.7, edgecolors='black', linewidth=0.3)
    
    ax_main.set_xscale('log')
    ax_main.set_yscale('log')
    ax_main.set_xlabel('Total chrR reads (chrR_A_n + chrR_B_n)', fontsize=12, fontweight='bold')
    ax_main.set_ylabel('Max chrR span (max of chrR_A_span, chrR_B_span)', fontsize=12, fontweight='bold')
    ax_main.grid(True, alpha=0.3, linestyle='--', which='both')
    
    # Top histogram (x-axis marginal) - log scale
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_top.hist(total_reads, bins=np.logspace(np.log10(min(total_reads)), 
                                                np.log10(max(total_reads)), 40),
                color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax_top.set_ylabel('Count', fontsize=10)
    ax_top.tick_params(labelbottom=False)
    ax_top.set_title('Ectopic rDNA Insertion Clusters (MAPQ > 20, Total reads ≥ 2)', 
                     fontsize=14, fontweight='bold', pad=10)
    
    # Right histogram (y-axis marginal)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
    # Use log-spaced bins for y-axis, handling zeros by adding small offset
    y_vals_nonzero = max_spans[max_spans > 0]
    if len(y_vals_nonzero) > 0:
        bins_y = np.logspace(np.log10(max(0.1, min(y_vals_nonzero))), 
                             np.log10(max(max_spans)), 40)
    else:
        bins_y = 40
    ax_right.hist(max_spans, bins=bins_y, orientation='horizontal',
                  color='coral', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax_right.set_xlabel('Count', fontsize=10)
    ax_right.tick_params(labelleft=False)
    
    # Add colorbar for density
    cbar_ax = fig.add_subplot(gs[1, 2])
    cbar = plt.colorbar(scatter, cax=cbar_ax)
    cbar.set_label('Point Density', fontsize=10)
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_hexbin_plot(total_reads, max_spans, output_path):
    """Create hexbin plot - great for showing density of overlapping points."""
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Filter out zero or negative spans for log scale
    mask = max_spans > 0
    filtered_reads = total_reads[mask]
    filtered_spans = max_spans[mask]
    
    if len(filtered_reads) == 0:
        print("Warning: No positive span values for hexbin plot")
        return
    
    hexbin = ax.hexbin(filtered_reads, filtered_spans, xscale='log', yscale='log',
                       gridsize=30, cmap='YlOrRd', mincnt=1, 
                       edgecolors='black', linewidths=0.2)
    
    ax.set_xlabel('Total chrR reads (chrR_A_n + chrR_B_n)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Max chrR span (max of chrR_A_span, chrR_B_span)', fontsize=12, fontweight='bold')
    ax.set_title('Ectopic rDNA Insertion Clusters - Hexbin Density Plot\n(MAPQ > 20, Total reads ≥ 2, Span > 0)', 
                 fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--', which='both')
    
    cb = plt.colorbar(hexbin, ax=ax)
    cb.set_label('Cluster Count per Bin', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def create_kde_plot(total_reads, max_spans, output_path):
    """Create 2D KDE contour plot."""
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Transform x to log scale for KDE
    log_reads = np.log10(total_reads)
    # Transform y to log scale for KDE (handle zeros)
    log_spans = np.log10(np.maximum(max_spans, 0.1))
    
    # Create KDE
    try:
        sns.kdeplot(x=log_reads, y=log_spans, cmap='viridis', fill=True, 
                    thresh=0, levels=20, alpha=0.6, ax=ax)
        
        # Overlay scatter points
        ax.scatter(log_reads, log_spans, s=30, alpha=0.5, 
                  edgecolors='black', linewidth=0.3, c='white')
        
        # Convert x-axis labels back to original scale
        ax.set_xlabel('Total chrR reads (chrR_A_n + chrR_B_n)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Max chrR span (max of chrR_A_span, chrR_B_span)', fontsize=12, fontweight='bold')
        ax.set_title('Ectopic rDNA Insertion Clusters - 2D Density Plot\n(MAPQ > 20, Total reads ≥ 2)', 
                     fontsize=14, fontweight='bold')
        
        # Set x-tick labels to show original values
        x_ticks = ax.get_xticks()
        ax.set_xticklabels([f'{10**x:.0f}' for x in x_ticks])
        
        # Set y-tick labels to show original values
        y_ticks = ax.get_yticks()
        ax.set_yticklabels([f'{10**y:.0f}' for y in y_ticks])
        
        ax.grid(True, alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"Could not create KDE plot: {e}")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze ectopic rDNA insertion clusters with enhanced visualizations'
    )
    parser.add_argument('directory', type=str,
                       help='Directory containing *_putative_bp.tsv files')
    parser.add_argument('-o', '--output', type=str, default='rDNA_clusters',
                       help='Output prefix for plots (default: rDNA_clusters)')
    parser.add_argument('--plot-type', type=str, default='all',
                       choices=['marginal', 'hexbin', 'kde', 'all'],
                       help='Type of plot to generate (default: all)')
    
    args = parser.parse_args()
    
    directory = Path(args.directory)
    if not directory.exists():
        print(f"Error: Directory '{directory}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    tsv_files = list(directory.glob('*_putative_bp.tsv'))
    
    if not tsv_files:
        print(f"Warning: No files matching '*_putative_bp.tsv' found in {directory}", 
              file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(tsv_files)} TSV file(s)")
    
    # Process all files
    all_data = []
    for filepath in sorted(tsv_files):
        print(f"Processing: {filepath.name}")
        data = process_file(filepath)
        all_data.extend(data)
        print(f"  - Found {len(data)} qualifying clusters")
    
    if not all_data:
        print("\nNo data points passed the filtering criteria!", file=sys.stderr)
        sys.exit(1)
    
    max_spans = np.array([d[0] for d in all_data])
    total_reads = np.array([d[1] for d in all_data])
    
    print(f"\nTotal qualifying clusters: {len(all_data)}")
    print(f"Max span range: {min(max_spans):.1f} - {max(max_spans):.1f}")
    print(f"Total reads range: {min(total_reads):.0f} - {max(total_reads):.0f}")
    print(f"Median reads: {np.median(total_reads):.0f}")
    print(f"Median span: {np.median(max_spans):.1f}")
    
    # Generate plots
    output_prefix = args.output
    
    if args.plot_type in ['marginal', 'all']:
        output_path = f"{output_prefix}_marginal.png"
        print(f"\nCreating marginal distribution plot...")
        create_marginal_plot(total_reads, max_spans, output_path)
        print(f"Saved: {output_path}")
    
    if args.plot_type in ['hexbin', 'all']:
        output_path = f"{output_prefix}_hexbin.png"
        print(f"Creating hexbin density plot...")
        create_hexbin_plot(total_reads, max_spans, output_path)
        print(f"Saved: {output_path}")
    
    if args.plot_type in ['kde', 'all']:
        output_path = f"{output_prefix}_kde.png"
        print(f"Creating 2D KDE plot...")
        create_kde_plot(total_reads, max_spans, output_path)
        print(f"Saved: {output_path}")
    
    # Save data to CSV
    csv_output = f"{output_prefix}_data.csv"
    df_output = pd.DataFrame(all_data, columns=['max_chrR_span', 'total_chrR_reads', 'source_file'])
    df_output.to_csv(csv_output, index=False)
    print(f"\nData saved to: {csv_output}")


if __name__ == '__main__':
    main()