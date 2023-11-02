configfile: "config/config.yaml"
SAMPLES = glob_wildcards(config["fastq"] + "{smp}.fq.gz").smp

# --- Build Rules --- #

## all                                : Build all final outputs
rule all:
    input:
        expand(config["results"] + "read_counts/totalcounts_mid/{smp}.totalcounts.tab",smp=SAMPLES),
        expand(config["results"] + "read_counts/totalcounts_full/{smp}.totalcounts.tab",smp=SAMPLES),
        expand(config["results"] + "read_counts/uniquecounts_mid/{smp}.uniquecounts.tab",smp=SAMPLES),
        expand(config["results"] + "read_counts/uniquecounts_full/{smp}.uniquecounts.tab",smp=SAMPLES),
        expand(config["results"] + "bigwig/{smp}.bw",smp=SAMPLES)


# --- Include Rules --- #
include: "rules/preprocess.smk"
include: "rules/align.smk"
include: "rules/map.smk"
include: "rules/bam2bw.smk"
include: "rules/tab.smk"
