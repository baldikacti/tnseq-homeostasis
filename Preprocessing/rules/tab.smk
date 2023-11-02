# --- Tab Files for Mid-Totalcounts --- #
## tab_generate_totalcounts           : Creates .tab files using trimmed totalcounts
rule tab_generate_totalcounts:
    input:
        script = config["script"] + "tab.R",
        data = config["results"] + "mapped/{smp}_middle_totalsum.txt"
    output:
        config["results"] + "read_counts/totalcounts_mid/{smp}.totalcounts.tab"
    shell:
        "Rscript {input.script} --data {input.data} --out {output}"

# --- Tab Files for Full-Totalcounts --- #
## tab_generate_totalcounts           : Creates .tab files using full length totalcounts
rule tab_generate_full_totalcounts:
    input:
        script = config["script"] + "tab.R",
        data = config["results"] + "mapped/{smp}_full_totalsum.txt"
    output:
        config["results"] + "read_counts/totalcounts_full/{smp}.totalcounts.tab"
    shell:
        "Rscript {input.script} --data {input.data} --out {output}"

# --- Tab Files for Mid-Uniquecounts --- #
## tab_generate_uniquecounts          : Creates .tab files using trimmed uniquecounts
rule tab_generate_uniquecounts:
    input:
        script = config["script"] + "tab.R",
        data = config["results"] + "mapped/{smp}_middle_uniquesum.txt"
    output:
        config["results"] + "read_counts/uniquecounts_mid/{smp}.uniquecounts.tab"
    shell:
        "Rscript {input.script} --data {input.data} --out {output}"

# --- Tab Files for Full-Uniquecounts --- #
## tab_generate_uniquecounts          : Creates .tab files using full length uniquecounts
rule tab_generate_full_uniquecounts:
    input:
        script = config["script"] + "tab.R",
        data = config["results"] + "mapped/{smp}_full_uniquesum.txt"
    output:
        config["results"] + "read_counts/uniquecounts_full/{smp}.uniquecounts.tab"
    shell:
        "Rscript {input.script} --data {input.data} --out {output}"