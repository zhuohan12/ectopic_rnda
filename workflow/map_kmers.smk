import os

output_dir = config["kmer_map"]["output_dir"]
# Make sure output_dir ends with slash (optional, for safety)
if not output_dir.endswith(os.sep):
    output_dir += os.sep


rule jellyfish_count:
    input:
        reads = config["kmer_map"]["jellyfish"]["rdna_path"]
    output:
        jf = os.path.join(config["kmer_map"]["output_dir"], f"mer_counts_k{config['kmer_map']['jellyfish']['k']}.jf")
    params:
        k = config["kmer_map"]["jellyfish"]["k"],
        s = config["kmer_map"]["jellyfish"]["hash_size"],
        t = config["kmer_map"]["jellyfish"]["threads"],
        outdir = os.path.dirname(
            os.path.join(config["kmer_map"]["output_dir"], f"mer_counts_k{config['kmer_map']['jellyfish']['k']}.jf")
        )
    shell:
        """
        mkdir -p {params.outdir}
        jellyfish count -m {params.k} -s {params.s} -t {params.t} -C {input.reads} -o {output.jf}
        """


rule jellyfish_dump:
    input:
        jf = os.path.join(output_dir, f"mer_counts_k{config['kmer_map']['jellyfish']['k']}.jf")
    output:
        fa = os.path.join(output_dir, f"mer_counts_k{config['kmer_map']['jellyfish']['k']}.fa")
    params:
        outdir = os.path.dirname(
            os.path.join(output_dir, f"mer_counts_k{config['kmer_map']['jellyfish']['k']}.fa")
        )
    shell:
        """
        mkdir -p {params.outdir}
        jellyfish dump {input.jf} > {output.fa}
        """

rule minimap:
    input:
        kmers = os.path.join(output_dir, f"mer_counts_k{config['kmer_map']['jellyfish']['k']}.fa"),
        ref = lambda wildcards: ref_dict[wildcards.ref]
    output:
        sam = os.path.join(output_dir, "{ref}_kmer_aligned.sam")
    params:
        outdir = os.path.dirname(
            os.path.join(output_dir, "{ref}_kmer_aligned.sam")
        )
    shell:
        """
        mkdir -p {params.outdir}
        minimap2 -a {input.ref} {input.kmers} > {output.sam}
        """

rule convert_bam:
    
