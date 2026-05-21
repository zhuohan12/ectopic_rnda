# Snakefile
configfile: "config.yaml"

# Include the main workflow files
include: "workflow/generate_kmer_set.smk"
include: "workflow/run_sr_pipeline_fastq.smk"

