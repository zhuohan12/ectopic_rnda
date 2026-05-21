#!/usr/bin/env python3
import sys, pysam
from collections import namedtuple

# ---- Tunables ----
L_MIN = 20       # minimum soft-clip length
MIN_MAPQ = 20    # minimum mapping quality of aligned part
MIN_MEAN_Q = 20  # mean Phred within the clip

ClipRec = namedtuple(
    "ClipRec",
    "qname chrom pos strand cigar mapq is_r1 side seq qual sa"
)

def phred_mean(qstr):
    return 0 if not qstr else sum(ord(c) - 33 for c in qstr) / len(qstr)

def softclip_ops(read):
    """Return [('L', len, start, end), ('R', len, start, end)] for left/right clips."""
    ct = read.cigartuples or []
    out = []
    n = len(ct)
    if n == 0:
        return out
    seqlen = read.query_length or 0
    if ct[0][0] == 4:  # left clip
        L = ct[0][1]
        out.append(("L", L, 0, L))
    if ct[-1][0] == 4:  # right clip
        R = ct[-1][1]
        out.append(("R", R, seqlen - R, seqlen))
    return out

def pick_single_clip(read, min_len=L_MIN):
    """Pick exactly one long soft clip (≥ min_len). Drop double-clipped (S..M..S) reads."""
    clips = softclip_ops(read)
    if len(clips) >= 2:
        return None
    long = [c for c in clips if c[1] >= min_len]
    if len(long) != 1:
        return None
    return long[0]  # (side, length, start, end)

def trim_lowq(seq, qual, qmin=MIN_MEAN_Q):
    """Drop clip if mean Q < qmin."""
    return (seq, qual) if phred_mean(qual) >= qmin else ("", "")

def is_valid(read):
    return (
        not read.is_secondary
        and not read.is_supplementary
        and read.mapping_quality >= MIN_MAPQ
    )

def clip_pos(read, side):
    """Coordinate at the clipped boundary."""
    return read.reference_start if side == "L" else read.reference_end

def revcomp(seq):
    """Reverse complement a DNA sequence."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def extract_clips(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    clips = []
    for r in bam.fetch(until_eof=True):
        if not is_valid(r):
            continue

        pick = pick_single_clip(r, L_MIN)
        if not pick:
            continue

        side, length, i0, i1 = pick
        seq = r.query_sequence[i0:i1]
        qual = (r.qual or "")[i0:i1]
        seq, qual = trim_lowq(seq, qual, MIN_MEAN_Q)
        if len(seq) < L_MIN:
            continue

        # filter out low-complexity clips (single base dominating)
        top = max(seq.count(b) for b in "ACGTNacgtn")
        if top / len(seq) > 0.8:
            continue

        # Reverse complement sequence if read maps to negative strand
        if r.is_reverse:
            seq = revcomp(seq)
            qual = qual[::-1]

        clips.append(ClipRec(
            qname=r.query_name,
            chrom=r.reference_name,
            pos=clip_pos(r, side),
            strand="-" if r.is_reverse else "+",
            cigar=r.cigarstring or "",
            mapq=r.mapping_quality,
            is_r1=r.is_read1,
            side=side,
            seq=seq,
            qual=qual,
            sa=r.get_tag("SA") if r.has_tag("SA") else ""
        ))
    bam.close()
    return clips

def write_fasta(clips, path):
    with open(path, "w") as fh:
        for c in clips:
            header = f">{c.qname}|{c.chrom}:{c.pos}|{c.side}|{c.strand}|r{int(c.is_r1)}"
            fh.write(header + "\n")
            fh.write(c.seq + "\n")

def write_tsv(clips, path):
    with open(path, "w") as fh:
        fh.write("qname\tchrom\tpos\tstrand\tside\tmapq\tcigar\tclip_len\n")
        for c in clips:
            fh.write(
                f"{c.qname}\t{c.chrom}\t{c.pos}\t{c.strand}\t{c.side}\t"
                f"{c.mapq}\t{c.cigar}\t{len(c.seq)}\n"
            )

def main(bam_in, fasta_out, tsv_out):
    clips = extract_clips(bam_in)
    write_fasta(clips, fasta_out)
    write_tsv(clips, tsv_out)
    sys.stderr.write(
        f"Extracted {len(clips)} single-clipped reads → {fasta_out}\n"
    )

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(
            "Usage: extract_soft_clips.py input.bam output.fasta output.tsv\n"
        )
    bam_in    = sys.argv[1]
    fasta_out = sys.argv[2]
    tsv_out   = sys.argv[3]
    main(bam_in, fasta_out, tsv_out)
