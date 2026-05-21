# Load configs
intermediate_dir  = config["dirs"]["intermediate"]
output_dir        = config["dirs"]["output"]
samples           = config["samples"]
input_ref         = config["sr_pipe"]["ref_path"]
kmer_list         = config["sr_pipe"]["kmer_list"]
sr_dir            = config["sr_pipe"]["sr_dir"]
min_kmers         = config["sr_pipe"]["min_kmers"]
import os
_scripts          = os.path.join(workflow.basedir, "..", "scripts", "main")
printf_exec       = os.path.join(_scripts, "01_parse_reads")
extract_cases     = os.path.join(_scripts, "02_extract_caseABC.py")
extract_soft_clips = os.path.join(_scripts, "03_extract_soft_clips.py")
find_candidates   = os.path.join(_scripts, "04_cluster_sr_insertions.py")
filter_candidates = os.path.join(_scripts, "05_filter_putative_bp.py")
num_threads       = config["sr_pipe"]["num_threads"]

wildcard_constraints:
    sample="SRR[0-9]+"

rule all_two:
    input:
        expand(f"{output_dir}/{{sample}}_putative_bp_filt.tsv", sample=samples),

#############################################
# Step 1: Scan reads using parse_reads
#############################################
rule scan_reads_r1:
    input:
        kmers=kmer_list,
        r1=lambda wc: f"{sr_dir}/{wc.sample}_1.fastq.gz"   # !! CHANGED: direct fastq.gz input
    output:
        tsv=f"{intermediate_dir}/readscan/{{sample}}_1_k21.tsv"
    params:
        exe=printf_exec,
        t=num_threads
    shell:
        """
        zcat {input.r1} | {params.exe} {input.kmers} /dev/stdin {output.tsv} {params.t}
        """

rule scan_reads_r2:
    input:
        kmers=f"{output_dir}/unique_rdna_k21.txt",
        r2=lambda wc: f"{sr_dir}/{wc.sample}_2.fastq.gz"   # !! CHANGED: direct fastq.gz input
    output:
        tsv=f"{intermediate_dir}/readscan/{{sample}}_2_k21.tsv"
    params:
        exe=printf_exec,
        t=num_threads
    shell:
        """
        zcat {input.r2} | {params.exe} {input.kmers} /dev/stdin {output.tsv} {params.t}
        """
        # !! CHANGED: zcat replaces the samtools fastq + named pipe approach

#############################################
# Step 2: Filter read names
#############################################
rule filter_read_kmers:
    input:
        r1=f"{intermediate_dir}/readscan/{{sample}}_1_k21.tsv",
        r2=f"{intermediate_dir}/readscan/{{sample}}_2_k21.tsv"
    output:
        reads=f"{intermediate_dir}/{{sample}}_k21_min{min_kmers}_readnames.txt"
    params:
        n=min_kmers,
        t=num_threads
    shell:
        r"""
        awk '$NF >= {params.n} {{print $1}}' {input.r1} {input.r2} | \
        sort --parallel={params.t} -u > {output.reads}
        """

#############################################
# Step 3: Create a subset of FASTQ file (BOTH R1 and R2 together)
#############################################
rule subset_fastq_both:
    input:
        kmer_list=f"{intermediate_dir}/{{sample}}_k21_min{min_kmers}_readnames.txt",
        r1=lambda wc: f"{sr_dir}/{wc.sample}_1.fastq.gz",  # !! CHANGED: direct fastq.gz input
        r2=lambda wc: f"{sr_dir}/{wc.sample}_2.fastq.gz"   # !! CHANGED: direct fastq.gz input
    output:
        r1=f"{intermediate_dir}/{{sample}}_1_filt.fastq",
        r2=f"{intermediate_dir}/{{sample}}_2_filt.fastq"
    params:
        t=num_threads
    shell:
        """
        seqtk subseq {input.r1} {input.kmer_list} > {output.r1} &
        seqtk subseq {input.r2} {input.kmer_list} > {output.r2} &
        wait
        """

#############################################
# Step 4: Mapping of subsetted FASTQ to CHM13
#############################################
rule map_fastqs:
    input:
        r1=f"{intermediate_dir}/{{sample}}_1_filt.fastq",
        r2=f"{intermediate_dir}/{{sample}}_2_filt.fastq",
        ref=input_ref
    output:
        bam=f"{intermediate_dir}/{{sample}}.bam"
    params:
        t=num_threads
    shell:
        "bwa-mem2 mem -t {params.t} {input.ref} {input.r1} {input.r2} | samtools view -bS - > {output.bam}"

#############################################
# Step 5: Filter out unmapped reads
#############################################
rule filter_fastq:
    input:
        bam_orig=f"{intermediate_dir}/{{sample}}.bam"
    output:
        bam_filt=f"{intermediate_dir}/{{sample}}_filt.bam"
    params:
        t=num_threads
    shell:
        "samtools view -@ {params.t} -F 3852 {input.bam_orig} -bo {output.bam_filt}"

#############################################
# Step 6: Extraction of cases A,B,C
#############################################
rule extract_cases:
    input:
        bam=f"{intermediate_dir}/{{sample}}_filt.bam"
    output:
        caseA=f"{intermediate_dir}/{{sample}}_caseA.bam",
        caseB=f"{intermediate_dir}/{{sample}}_caseB.bam",
        caseC=f"{intermediate_dir}/{{sample}}_caseC.bam"
    params:
        script=extract_cases
    shell:
        "python3 {params.script} {input.bam} {output.caseA} {output.caseB} {output.caseC}"

#############################################
# Step 7: Extraction of soft clips
#############################################
rule extract_cases_A:
    input:
        caseA=f"{intermediate_dir}/{{sample}}_caseA.bam"
    output:
        fasta=f"{intermediate_dir}/{{sample}}_soft_clips_A.fasta",
        tsv=f"{intermediate_dir}/{{sample}}_soft_clips_A.tsv"
    params:
        script=extract_soft_clips
    shell:
        "python3 {params.script} {input.caseA} {output.fasta} {output.tsv}"

rule extract_cases_B:
    input:
        caseB=f"{intermediate_dir}/{{sample}}_caseB.bam"
    output:
        fasta=f"{intermediate_dir}/{{sample}}_soft_clips_B.fasta",
        tsv=f"{intermediate_dir}/{{sample}}_soft_clips_B.tsv"
    params:
        script=extract_soft_clips
    shell:
        "python3 {params.script} {input.caseB} {output.fasta} {output.tsv}"

#############################################
# Step 8: Remap soft clips
#############################################
rule remap_soft_clips_A:
    input:
        fasta_A=f"{intermediate_dir}/{{sample}}_soft_clips_A.fasta",
        ref=input_ref
    output:
        bam_A=f"{intermediate_dir}/{{sample}}_soft_clips_A_mapped.bam"
    params:
        t=num_threads
    shell:
        "bwa-mem2 mem -t {params.t} -k 11 -T 10 -a {input.ref} {input.fasta_A} | samtools view -bS - > {output.bam_A}"

rule remap_soft_clips_B:
    input:
        fasta_B=f"{intermediate_dir}/{{sample}}_soft_clips_B.fasta",
        ref=input_ref
    output:
        bam_B=f"{intermediate_dir}/{{sample}}_soft_clips_B_mapped.bam"
    params:
        t=num_threads
    shell:
        "bwa-mem2 mem -t {params.t} -k 11 -T 10 -a {input.ref} {input.fasta_B} | samtools view -bS - > {output.bam_B}"

#############################################
# Step 9: Cluster breakpoints
#############################################
rule cluster_breakpoints:
    input:
        bam_A=f"{intermediate_dir}/{{sample}}_soft_clips_A_mapped.bam",
        bam_B=f"{intermediate_dir}/{{sample}}_soft_clips_B_mapped.bam",
        bam_C=f"{intermediate_dir}/{{sample}}_caseC.bam"
    output:
        tsv=f"{intermediate_dir}/{{sample}}_putative_bp.tsv"
    params:
        script=find_candidates
    shell:
        "python3 {params.script} {input.bam_A} {input.bam_B} {input.bam_C} -o {output.tsv}"

#############################################
# Step 10: Filter breakpoints
#############################################
rule filt_breakpoints:
    input:
        all_bp=f"{intermediate_dir}/{{sample}}_putative_bp.tsv",
    output:
        tsv=f"{output_dir}/{{sample}}_putative_bp_filt.tsv"
    params:
        script=filter_candidates
    shell:
        "python3 {params.script} {input.all_bp} {output.tsv}"
