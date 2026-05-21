#!/usr/bin/env python3
"""
cluster_rdna_insertions.py

Goal: Cluster ectopic rDNA insertion signals from A/B/C read sets.
      Output clusters with per-case breakpoint coordinates.

Key features:
  - Outputs per-case genome coordinates: A_min, A_max, A_span, B_min, B_max, B_span, etc.
  - Outputs per-case rDNA coordinates: chrR_A_min, chrR_A_max, chrR_A_span, etc.
  - Outputs chrR_C_mad: median absolute deviation of Case C rDNA positions (robust span metric)
    MAD is only computed when chrR_C_n >= 3; set to '.' otherwise (insufficient data).
  - For Case C: extracts chrR position from discordant mate mapped to rDNA (FIXED)
  - Minimal pre-filtering: total >= min_total (default: 3) to remove obvious noise

PRE-FILTERS (technical, not biological):
  - min_mapq >= 0 (default: 0, can be raised to filter multi-mappers)
  - total >= 3 (default: 3, removes 1-2 read technical noise)

These are NOT quality filters - they remove obvious technical artifacts.
Quality filtering (nB, B/A ratio, etc.) should be done downstream.

Strand-aware breakpoint coordinates (junction-facing edge):
  - Case A and C: forward -> reference_end, reverse -> reference_start
  - Case B:       forward -> reference_start, reverse -> reference_end

Input: three BAMs (A/B/C read sets)
Output: TSV of merged loci with per-case statistics
"""

import sys
import argparse
import pysam
from collections import namedtuple

Hit = namedtuple("Hit", "chrom pos strand mapq rdna_pos")

def parse_chrR_pos(qname: str):
    """
    Parse chrR coordinate from read name token like: ...|chrR:13573|...
    Returns int or None if not parseable.
    """
    try:
        for part in qname.split("|"):
            if part.startswith("chrR:"):
                return int(part.split(":", 1)[1])
    except Exception:
        pass
    return None

def pick_breakpoint(aln, case: str):
    """
    Strand-aware breakpoint selection (junction-facing edge).

    Case A and C:
      + strand -> reference_end
      - strand -> reference_start

    Case B:
      + strand -> reference_start
      - strand -> reference_end
    """
    if case in ("A", "C"):
        return aln.reference_start if aln.is_reverse else aln.reference_end
    if case == "B":
        return aln.reference_end if aln.is_reverse else aln.reference_start
    raise ValueError(f"Unknown case: {case}")

def compute_mad(values):
    """
    Compute median absolute deviation (MAD) from a sorted list of numeric values.
    MAD = median(|xi - median(x)|)
    Returns None if fewer than 3 values (insufficient data for robust estimate).
    """
    n = len(values)
    if n < 3:
        return None
    median = values[n // 2]
    abs_devs = sorted(abs(v - median) for v in values)
    return abs_devs[n // 2]

def collect_hits(bam_path: str, case: str, min_mapq: int):
    """
    Collect hits from BAM file.

    For Case A/B: extract chrR position from read name (soft-clip remapped)
    For Case C: extract chrR position from MATE alignment (discordant pair)
    """
    hits = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue
            if aln.mapping_quality < min_mapq:
                continue
            if aln.is_secondary:
                continue

            chrom = aln.reference_name
            strand = "-" if aln.is_reverse else "+"
            pos = pick_breakpoint(aln, case)

            if pos is None:
                continue

            # ============ FIXED SECTION FOR CASE C ============
            rdna_pos = None

            if case in ("A", "B"):
                # Cases A/B: chrR position from read name (soft-clip remapped)
                rdna_pos = parse_chrR_pos(aln.query_name)

            elif case == "C":
                # Case C: chrR position from MATE alignment
                if not aln.mate_is_unmapped:
                    # Get mate's reference name
                    mate_chrom = aln.next_reference_name

                    # Check if mate maps to chrR
                    if mate_chrom == "chrR":
                        # Extract mate's position on chrR
                        rdna_pos = aln.next_reference_start
            # ==================================================

            hits.append(Hit(chrom, pos, strand, aln.mapping_quality, rdna_pos))
    return hits

def cluster_hits(hits, cluster_dist: int):
    """
    Cluster hits within cluster_dist per (chrom, strand).
    Returns list of cluster dicts with per-cluster statistics.
    """
    if not hits:
        return []

    hits.sort(key=lambda h: (h.chrom, h.strand, h.pos))
    clusters = []
    cur = [hits[0]]

    for h in hits[1:]:
        last = cur[-1]
        if h.chrom == last.chrom and h.strand == last.strand and (h.pos - last.pos) <= cluster_dist:
            cur.append(h)
        else:
            clusters.append(cur)
            cur = [h]
    clusters.append(cur)

    out = []
    for c in clusters:
        chrom = c[0].chrom
        strand = c[0].strand

        # Genome-side coordinates
        positions = sorted(h.pos for h in c)
        n = len(c)
        pos_min = positions[0]
        pos_max = positions[-1]
        pos_span = pos_max - pos_min
        pos_med = positions[n // 2]

        # MAPQ
        mean_mapq = sum(h.mapq for h in c) / n

        # rDNA-side coordinates (chrR)
        rdna_vals = sorted([h.rdna_pos for h in c if h.rdna_pos is not None])
        if rdna_vals:
            rdna_n = len(rdna_vals)
            rdna_min = rdna_vals[0]
            rdna_max = rdna_vals[-1]
            rdna_span = rdna_max - rdna_min
            rdna_med = rdna_vals[rdna_n // 2]
        else:
            rdna_n = 0
            rdna_min = None
            rdna_max = None
            rdna_span = None
            rdna_med = None

        out.append(dict(
            chrom=chrom,
            strand=strand,
            pos_min=pos_min,
            pos_max=pos_max,
            pos_span=pos_span,
            pos_med=pos_med,
            n=n,
            mean_mapq=mean_mapq,
            rdna_n=rdna_n,
            rdna_min=rdna_min,
            rdna_max=rdna_max,
            rdna_span=rdna_span,
            rdna_med=rdna_med,
            hits=c,  # Keep raw hits for exact locus-level stats
        ))
    return out

def merge_clusters_across_cases(Ac, Bc, Cc, merge_dist: int, min_total: int):
    """
    Merge clusters across A/B/C cases into loci based on genomic proximity.

    Output per-case coordinates:
      - Genome: A_min, A_max, A_span, B_min, B_max, B_span, C_min, C_max, C_span
      - rDNA: chrR_A_min, chrR_A_max, chrR_A_span, chrR_B_min, ..., chrR_C_span
      - chrR_C_mad: MAD of Case C rDNA positions (computed when chrR_C_n >= 3)

    Pre-filter: total >= min_total (removes 1-2 read noise)
    """
    allc = []
    for case, clusters in (("A", Ac), ("B", Bc), ("C", Cc)):
        for cl in clusters:
            allc.append((case, cl))

    # Group by chromosome
    by_chrom = {}
    for case, cl in allc:
        by_chrom.setdefault(cl["chrom"], []).append((case, cl))

    merged = []

    for chrom, items in by_chrom.items():
        # Sort by start position
        items.sort(key=lambda x: x[1]["pos_min"])

        cur = [items[0]]
        cur_start = items[0][1]["pos_min"]
        cur_end = items[0][1]["pos_max"]

        def flush(group, gstart, gend):
            """Process a merged locus from grouped clusters"""

            # Separate clusters by case
            clusters_A = [cl for case, cl in group if case == "A"]
            clusters_B = [cl for case, cl in group if case == "B"]
            clusters_C = [cl for case, cl in group if case == "C"]

            # Per-case read counts
            nA = sum(cl["n"] for cl in clusters_A)
            nB = sum(cl["n"] for cl in clusters_B)
            nC = sum(cl["n"] for cl in clusters_C)
            total = nA + nB + nC

            # PRE-FILTER: Remove obvious noise (total < min_total)
            if total < min_total:
                return

            # Per-case genome coordinates
            def get_coords(clusters, default="."):
                if not clusters:
                    return default, default, default, default
                all_pos = []
                for cl in clusters:
                    all_pos.extend([h.pos for h in cl["hits"]])
                all_pos.sort()
                return (min(all_pos), max(all_pos), 
                       max(all_pos) - min(all_pos),
                       all_pos[len(all_pos)//2])

            A_min, A_max, A_span, A_med = get_coords(clusters_A)
            B_min, B_max, B_span, B_med = get_coords(clusters_B)
            C_min, C_max, C_span, C_med = get_coords(clusters_C)

            # Per-case rDNA coordinates
            def get_rdna_coords(clusters, default="."):
                if not clusters:
                    return default, default, default, default, 0, default
                all_rdna = []
                for cl in clusters:
                    all_rdna.extend([h.rdna_pos for h in cl["hits"] if h.rdna_pos is not None])
                if not all_rdna:
                    return default, default, default, default, 0, default
                all_rdna.sort()
                rdna_mad = compute_mad(all_rdna)  # None if n < 3
                rdna_mad_out = rdna_mad if rdna_mad is not None else default
                return (min(all_rdna), max(all_rdna),
                       max(all_rdna) - min(all_rdna),
                       all_rdna[len(all_rdna)//2],
                       len(all_rdna),
                       rdna_mad_out)

            chrR_A_min, chrR_A_max, chrR_A_span, chrR_A_med, chrR_A_n, chrR_A_mad = get_rdna_coords(clusters_A)
            chrR_B_min, chrR_B_max, chrR_B_span, chrR_B_med, chrR_B_n, chrR_B_mad = get_rdna_coords(clusters_B)
            chrR_C_min, chrR_C_max, chrR_C_span, chrR_C_med, chrR_C_n, chrR_C_mad = get_rdna_coords(clusters_C)

            # Modal strand per case
            def get_strand(clusters):
                if not clusters:
                    return "."
                from collections import Counter
                strands = [cl["strand"] for cl in clusters for _ in range(cl["n"])]
                if not strands:
                    return "."
                return Counter(strands).most_common(1)[0][0]

            strand_A = get_strand(clusters_A)
            strand_B = get_strand(clusters_B)
            strand_C = get_strand(clusters_C)

            # Weighted mean MAPQ
            all_hits = []
            for _, cl in group:
                all_hits.extend(cl["hits"])
            mean_mapq = sum(h.mapq for h in all_hits) / len(all_hits) if all_hits else 0.0

            merged.append(dict(
                chrom=chrom,
                start=gstart,
                end=gend,
                span=gend - gstart,

                # Per-case read counts
                nA=nA,
                nB=nB,
                nC=nC,
                total=total,

                # Per-case strands
                strand_A=strand_A,
                strand_B=strand_B,
                strand_C=strand_C,

                # Per-case genome coordinates
                A_min=A_min, A_max=A_max, A_span=A_span, A_med=A_med,
                B_min=B_min, B_max=B_max, B_span=B_span, B_med=B_med,
                C_min=C_min, C_max=C_max, C_span=C_span, C_med=C_med,

                # Per-case rDNA coordinates
                chrR_A_n=chrR_A_n,
                chrR_A_min=chrR_A_min, chrR_A_max=chrR_A_max,
                chrR_A_span=chrR_A_span, chrR_A_med=chrR_A_med,
                chrR_A_mad=chrR_A_mad,

                chrR_B_n=chrR_B_n,
                chrR_B_min=chrR_B_min, chrR_B_max=chrR_B_max,
                chrR_B_span=chrR_B_span, chrR_B_med=chrR_B_med,
                chrR_B_mad=chrR_B_mad,

                chrR_C_n=chrR_C_n,
                chrR_C_min=chrR_C_min, chrR_C_max=chrR_C_max,
                chrR_C_span=chrR_C_span, chrR_C_med=chrR_C_med,
                chrR_C_mad=chrR_C_mad,

                # Overall MAPQ
                mean_mapq=mean_mapq,
            ))

        # Merge clusters within merge_dist
        for case, cl in items[1:]:
            if cl["pos_min"] <= cur_end + merge_dist:
                cur.append((case, cl))
                if cl["pos_max"] > cur_end:
                    cur_end = cl["pos_max"]
            else:
                flush(cur, cur_start, cur_end)
                cur = [(case, cl)]
                cur_start = cl["pos_min"]
                cur_end = cl["pos_max"]

        flush(cur, cur_start, cur_end)

    return merged

def main():
    ap = argparse.ArgumentParser(
        description="Cluster rDNA ectopic insertion signals with minimal pre-filtering."
    )
    ap.add_argument("bamA", help="Case A BAM")
    ap.add_argument("bamB", help="Case B BAM")
    ap.add_argument("bamC", help="Case C BAM")
    ap.add_argument("-o", "--out", default="clusters_all.tsv", help="Output TSV")

    # Technical pre-filters (not quality filters)
    ap.add_argument("--min-mapq", type=int, default=0, 
                    help="Minimum MAPQ for read collection (default: 0)")
    ap.add_argument("--min-total", type=int, default=3,
                    help="Minimum total reads per locus (default: 3)")
    # Clustering parameters
    ap.add_argument("--cluster-dist", type=int, default=50, 
                    help="Clustering distance within case (bp, default: 50)")
    ap.add_argument("--merge-dist", type=int, default=100, 
                    help="Distance to merge clusters across cases (bp, default: 100)")

    args = ap.parse_args()

    # Collect hits
    sys.stderr.write("Collecting hits from BAMs...\n")
    hitsA = collect_hits(args.bamA, case="A", min_mapq=args.min_mapq)
    hitsB = collect_hits(args.bamB, case="B", min_mapq=args.min_mapq)
    hitsC = collect_hits(args.bamC, case="C", min_mapq=args.min_mapq)

    sys.stderr.write(f"  Case A: {len(hitsA)} hits\n")
    sys.stderr.write(f"  Case B: {len(hitsB)} hits\n")
    sys.stderr.write(f"  Case C: {len(hitsC)} hits\n")

    # Cluster within each case
    sys.stderr.write("Clustering within cases...\n")
    clustersA = cluster_hits(hitsA, args.cluster_dist)
    clustersB = cluster_hits(hitsB, args.cluster_dist)
    clustersC = cluster_hits(hitsC, args.cluster_dist)

    sys.stderr.write(f"  Case A: {len(clustersA)} clusters\n")
    sys.stderr.write(f"  Case B: {len(clustersB)} clusters\n")
    sys.stderr.write(f"  Case C: {len(clustersC)} clusters\n")

    # Merge across cases
    sys.stderr.write(f"Merging clusters across cases (min_total >= {args.min_total})...\n")
    merged = merge_clusters_across_cases(
        clustersA, clustersB, clustersC,
        merge_dist=args.merge_dist,
        min_total=args.min_total
    )

    sys.stderr.write(f"Merged into {len(merged)} loci (total >= {args.min_total})\n")

    # Write output
    with open(args.out, "w") as out:
        # Header
        header = [
            "chrom", "start", "end", "span",
            "nA", "nB", "nC", "total",
            "strand_A", "strand_B", "strand_C",
            "A_min", "A_max", "A_span", "A_med",
            "B_min", "B_max", "B_span", "B_med",
            "C_min", "C_max", "C_span", "C_med",
            "chrR_A_n", "chrR_A_min", "chrR_A_max", "chrR_A_span", "chrR_A_med", "chrR_A_mad",
            "chrR_B_n", "chrR_B_min", "chrR_B_max", "chrR_B_span", "chrR_B_med", "chrR_B_mad",
            "chrR_C_n", "chrR_C_min", "chrR_C_max", "chrR_C_span", "chrR_C_med", "chrR_C_mad",
            "mean_mapq"
        ]
        out.write("\t".join(header) + "\n")

        # Data rows
        for m in merged:
            row = [
                m["chrom"], m["start"], m["end"], m["span"],
                m["nA"], m["nB"], m["nC"], m["total"],
                m["strand_A"], m["strand_B"], m["strand_C"],
                m["A_min"], m["A_max"], m["A_span"], m["A_med"],
                m["B_min"], m["B_max"], m["B_span"], m["B_med"],
                m["C_min"], m["C_max"], m["C_span"], m["C_med"],
                m["chrR_A_n"], m["chrR_A_min"], m["chrR_A_max"], m["chrR_A_span"], m["chrR_A_med"], m["chrR_A_mad"],
                m["chrR_B_n"], m["chrR_B_min"], m["chrR_B_max"], m["chrR_B_span"], m["chrR_B_med"], m["chrR_B_mad"],
                m["chrR_C_n"], m["chrR_C_min"], m["chrR_C_max"], m["chrR_C_span"], m["chrR_C_med"], m["chrR_C_mad"],
                f"{m['mean_mapq']:.1f}"
            ]
            out.write("\t".join(str(x) for x in row) + "\n")

    sys.stderr.write(f"Wrote {args.out}\n")

if __name__ == "__main__":
    main()