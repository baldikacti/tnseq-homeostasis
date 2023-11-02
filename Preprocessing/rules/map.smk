# --- Assign Genome Positions --- #
## genomecov                          : Assigns genome positions to counts
rule genomecov:
    input:
        config["results"] + "bwa_aln/{smp}_marked.bam"
    output:
        config["results"] + "mapped/{smp}_count.txt"
    shell:
        "bedtools genomecov -5 -bg -ibam {input} > {output}"

# --- Map to Mid-Totalcounts --- #
## mapping_totalcounts                : Creates count files with trimmed totalcounts
rule mapping_totalcounts:
    input:
        ref= config["ref"]["map_mid"],
        ct= config["results"] + "mapped/{smp}_count.txt"
    output:
        config["results"] + "mapped/{smp}_middle_totalsum.txt"
    shell:
        "bedtools map -a {input.ref} -b {input.ct} -c 4 -o sum > {output}"

# --- Map to Full-Totalcounts --- #
## mapping_totalcounts                : Creates count files with full length totalcounts
rule mapping_full_totalcounts:
    input:
        ref= config["ref"]["map_full"],
        ct= config["results"] + "mapped/{smp}_count.txt"
    output:
        config["results"] + "mapped/{smp}_full_totalsum.txt"
    shell:
        "bedtools map -a {input.ref} -b {input.ct} -c 4 -o sum > {output}"

# --- Map to Mid-Uniquecounts --- #
## mapping_uniquecounts               : Creates count files with trimmed uniquecounts
rule mapping_uniquecounts:
    input:
        ref= config["ref"]["map_mid"],
        ct= config["results"] + "mapped/{smp}_count.txt"
    output:
        config["results"] + "mapped/{smp}_middle_uniquesum.txt"
    shell:
        "bedtools map -a {input.ref} -b {input.ct} -c 4 -o count > {output}"

# --- Map to Full-Uniquecounts --- #
## mapping_uniquecounts               : Creates count files with full length uniquecounts
rule mapping_full_uniquecounts:
    input:
        ref= config["ref"]["map_full"],
        ct= config["results"] + "mapped/{smp}_count.txt"
    output:
        config["results"] + "mapped/{smp}_full_uniquesum.txt"
    shell:
        "bedtools map -a {input.ref} -b {input.ct} -c 4 -o count > {output}"