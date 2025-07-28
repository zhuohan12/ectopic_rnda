# Snakefile
configfile: "/d/hd09/zhuohan/ectopic_rdna_proj/config/config.yaml"

# Include the main workflow files
include: "/d/hd09/zhuohan/ectopic_rdna_proj/workflow/map_kmers.smk"
# include: "/d/hd09/zhuohan/ectopic_rdna_proj/workflow/kat_comp.smk"

ref_names = config["kmer_map"]["minimap"]["ref_names"]
ref_paths = config["kmer_map"]["minimap"]["ref_paths"]
ref_dict = dict(zip(ref_names, ref_paths))

output_dir = config["kmer_map"]["output_dir"]
rule all:
    input:
        expand(os.path.join(output_dir, "{ref}_kmer_aligned.sam"), ref=ref_names)
