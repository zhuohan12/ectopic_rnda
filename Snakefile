# Snakefile
configfile: "/d/hd09/zhuohan/ectopic_rdna_proj/config.yaml"

# Include the main workflow files
include: "/d/hd09/zhuohan/ectopic_rdna_proj/workflow/generate_kmer_set.smk"
include: "/d/hd09/zhuohan/ectopic_rdna_proj/workflow/run_sr_pipeline.smk"

