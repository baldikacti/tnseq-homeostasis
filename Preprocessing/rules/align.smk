# --- BWA Alignment --- #
## bwa_mem                            : Creates .bam files from .fq.gz
rule bwa_mem:
    input:
        fa= config["ref"]["fasta"],
        fq= config["results"] + "preprocess/{smp}_pruned.fq.gz"
    output:
        config["results"] + "bwa_aln/{smp}.bam"
    threads: 4
    shell:
        "bwa mem -t {threads} {input.fa} {input.fq} | samtools sort -@2 - > {output}"

# --- Remove PCR Duplicates --- #
## rm_dupe                            : Removes PCR duplicates using JE
rule rm_dupe:
    input:
        bam = config["results"] + "bwa_aln/{smp}.bam",
        je = config["script"] + "je_1.2_bundle.jar"
    output:
        bam = config["results"] + "bwa_aln/{smp}_marked.bam",
        metrics = config["results"] + "bwa_aln/dedup_metrics/{smp}_metrics.txt"
    shell:
        "java -jar {input.je} markdupes I={input.bam} O={output.bam} M={output.metrics} MM=0 REMOVE_DUPLICATES=TRUE"