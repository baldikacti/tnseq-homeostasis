# --- Generate BigWig files --- #
## bam2bw                          : Creates BigWig files from deduplicated bam files
rule bam2bw:
    input:
        config["results"] + "bwa_aln/{smp}_marked.bam"
    output:
        bw = config["results"] + "bigwig/{smp}.bw",
        bai = config["results"] + "bwa_aln/{smp}_marked.bam.bai"
    conda:
        "envs/environment.yaml"
    threads: 2
    shell:
        """
        samtools index -@ {threads} {input}
        bamCoverage -b {input} -o {output.bw}
        """