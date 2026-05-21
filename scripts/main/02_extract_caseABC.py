#!/usr/bin/env python3
import pysam
from collections import defaultdict

RDNA_CONTIG = "chrR"

def clip_length(read):
    if read.cigartuples is None:
        return 0
    return sum(length for op, length in read.cigartuples if op == 4)

def has_sa_cross(read):
    if not read.has_tag("SA"):
        return False
    sa = read.get_tag("SA")
    for aln in sa.split(";"):
        if aln:
            contig = aln.split(",")[0]
            if (read.reference_name == RDNA_CONTIG) != (contig == RDNA_CONTIG):
                return True
    return False

def is_clean_anchor(read):
    if read.is_secondary or read.is_supplementary:
        return False
    if "S" in read.cigarstring or "H" in read.cigarstring:
        return False
    return True

def is_spanning(read, min_clip=20):
    if read.is_secondary:
        return False
    if clip_length(read) >= min_clip:
        return True
    if has_sa_cross(read):
        return True
    return False

def classify_pair(reads):
    """Return 'A', 'B', or 'C' for the pair type, or None."""
    primaries = [r for r in reads if not r.is_secondary]
    if len(primaries) < 2:
        return None

    r1, r2 = primaries[:2]
    r1_in_rdna = r1.reference_name == RDNA_CONTIG
    r2_in_rdna = r2.reference_name == RDNA_CONTIG

    # --- Case A ---
    if (
        (r1_in_rdna and is_clean_anchor(r1) and r2_in_rdna and is_spanning(r2))
        or (r2_in_rdna and is_clean_anchor(r2) and r1_in_rdna and is_spanning(r1))
    ):
        return "A"

    # --- Case B ---
    if (
        (r1_in_rdna and is_spanning(r1) and is_clean_anchor(r2) and not r2_in_rdna)
        or (r2_in_rdna and is_spanning(r2) and is_clean_anchor(r1) and not r1_in_rdna)
    ):
        return "B"

    # --- Case C ---
    if r1_in_rdna != r2_in_rdna and is_clean_anchor(r1) and is_clean_anchor(r2):
        return "C"

    return None


def main(in_bam, outA, outB, outC, report_every=500000, nthreads=16):
    bam_in = pysam.AlignmentFile(in_bam, "rb", threads=nthreads)
    bamA   = pysam.AlignmentFile(outA, "wb", template=bam_in, threads=nthreads)
    bamB   = pysam.AlignmentFile(outB, "wb", template=bam_in, threads=nthreads)
    bamC   = pysam.AlignmentFile(outC, "wb", template=bam_in, threads=nthreads)

    counts = {"A": 0, "B": 0, "C": 0}
    total = 0
    prev_qname, bucket = None, []

    for read in bam_in.fetch(until_eof=True):
        if prev_qname is not None and read.query_name != prev_qname:
            case = classify_pair(bucket)
            if case:
                for r in bucket:
                    if case == "A":
                        bamA.write(r)
                    elif case == "B":
                        bamB.write(r)
                    elif case == "C":
                        bamC.write(r)
                counts[case] += 1

            total += 1
            if total % report_every == 0:
                print(f"Processed {total:,} pairs")

            bucket = []

        bucket.append(read)
        prev_qname = read.query_name

    # finalize last pair
    if bucket:
        case = classify_pair(bucket)
        if case:
            for r in bucket:
                if case == "A":
                    bamA.write(r)
                elif case == "B":
                    bamB.write(r)
                elif case == "C":
                    bamC.write(r)
            counts[case] += 1
        total += 1

    bam_in.close()
    bamA.close()
    bamB.close()
    bamC.close()

    print("\nFinished.")
    for k in ["A", "B", "C"]:
        print(f"Wrote {counts[k]} case {k} pairs to case{k}.bam")
    print("Unclassified pairs:", total - sum(counts.values()))
    print("Total pairs processed:", total)


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        sys.exit("Usage: extract_cases.py input.bam [caseA.bam caseB.bam caseC.bam]")
    in_bam = sys.argv[1]
    outA = sys.argv[2] if len(sys.argv) > 2 else "caseA.bam"
    outB = sys.argv[3] if len(sys.argv) > 3 else "caseB.bam"
    outC = sys.argv[4] if len(sys.argv) > 4 else "caseC.bam"
    main(in_bam, outA, outB, outC)
