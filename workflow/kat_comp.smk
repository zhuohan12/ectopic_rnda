import os

k_values = config["kat"]["k_values"]
ref_path = config["kat"]["ref"]
rdna_path = config["kat"]["rdna"]
output_dir = config["kat"]["output_dir"]
hash_size = config["kat"]["hash_size"]
threads = config["kat"]["threads"]

# Clean reference genome name
ref_name = os.path.basename(ref_path).lower()
for ext in [".fasta", ".fa", ".fna"]:
    if ref_name.endswith(ext):
        ref_name = ref_name[:-len(ext)]
        break

rule kat_comp:
    input:
        ref=ref_path,
        rdna=rdna_path
    output:
        matrix=os.path.join(output_dir, ref_name, "kat.k{km}.comp.mx")
    params:
        outprefix=lambda wc: os.path.join(output_dir, ref_name, f"kat.k{wc.km}"),
        outdir=lambda wc: os.path.dirname(os.path.join(output_dir, ref_name, f"kat.k{wc.km}.comp.mx"))
    threads: threads
    log:
        os.path.join(output_dir, ref_name, "kat.k{km}.log")
    wildcard_constraints:
        km="\\d+"
    shell:
        """
        mkdir -p {params.outdir}
        kat comp -o {params.outprefix} -m {wildcards.km} -H {hash_size} -t {threads} {input.ref} {input.rdna} > {log} 2>&1
        """