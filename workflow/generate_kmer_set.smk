# Snakemake: rDNA‑unique k-mer mapper

#############################################
# Load config values
#############################################
config_k          = config["kmer_set"]["k"]
input_rdna        = config["kmer_set"]["rdna_path"]
input_ref         = config["kmer_set"]["ref_path"]
num_threads       = config["kmer_set"]["num_threads"]
output_dir        = config["dirs"]["output"]
intermediate_dir  = config["dirs"]["intermediate"]
output_dir        = config["dirs"]["output"]


#############################################
# Final target
#############################################
rule all:
    input:
        f"{output_dir}/unique_rdna_k{config_k}.txt"

#############################################
# Step 1: Count genome k-mers
#############################################
rule count_genome_kmers:
    input:
        fasta=input_ref
    output:
        jf=temp(f"{intermediate_dir}/genome_k{config_k}.jf")
    params:
        k=config_k,
        t=num_threads
    shell:
        """
        jellyfish count -m {params.k} -s 200M -t {params.t} \
            -C {input.fasta} -o {output.jf}
        """

#############################################
# Step 2: Count rDNA k-mers
#############################################
rule count_rdna_kmers:
    input:
        fasta=input_rdna
    output:
        jf=temp(f"{intermediate_dir}/rdna_k{config_k}.jf")
    params:
        k=config_k,
        t=num_threads
    shell:
        """
        jellyfish count -m {params.k} -s 200M -t {params.t} \
            -C {input.fasta} -o {output.jf}
        """

#############################################
# Step 3: Dump Jellyfish k-mers
#############################################
rule dump_genome_kmers:
    input:
        jf=f"{intermediate_dir}/genome_k{config_k}.jf"
    output:
        fa=f"{intermediate_dir}/genome_k{config_k}.fa"
    shell:
        "jellyfish dump -c {input.jf} > {output.fa}"

rule dump_rdna_kmers:
    input:
        jf=f"{intermediate_dir}/rdna_k{config_k}.jf"
    output:
        fa=f"{intermediate_dir}/rdna_k{config_k}.fa"
    shell:
        "jellyfish dump -c {input.jf} > {output.fa}"

#############################################
# Step 4: Subtract kmers (unique to rDNA)
#############################################
rule subtract_kmers:
    input:
        genome=f"{intermediate_dir}/genome_k{config_k}.fa",
        rdna=f"{intermediate_dir}/rdna_k{config_k}.fa"
    output:
        unique=f"{output_dir}/unique_rdna_k{config_k}.txt"
    shell:
        r"""
        awk 'NR==FNR {{seen[$1]=1; next}} !($1 in seen) {{print $1}}' \
            {input.genome} {input.rdna} > {output.unique}
        """
