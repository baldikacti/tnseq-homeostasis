# Chienlab-TNseq

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.32.4-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

## Introduction

**chienlab-tnseq** is a Snakemake pipeline for performing bacterial TNseq alignment and annotation.

This repository contains version of the pipeline used in the publication.

The maintained and uptodate version of the pipeline can be found [here](https://github.com/baldikacti/chienlab-tnseq).

## Pipeline summary

The pipeline will perform the following steps:

1. Removes reads that does not contain the static transposon region. ([`Seqkit`](https://bioinf.shenwei.me/seqkit/))
2. Removes PCR duplicates using UMIs. ([`JE`](https://github.com/gbcs-embl/Je))
3. Align reads to reference genome. ([`BWA-MEM`](https://github.com/lh3/bwa/))
4. Generates BigWig files for visualization in genome browsers ([`deeptools`](https://deeptools.readthedocs.io/en/develop/))
5. Assigns 5' position to counts. ([`bedtools genomocov`](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html))
6. Maps position counts to gene features using either clipped or full sized genes. ([`bedtools map`](https://bedtools.readthedocs.io/en/latest/content/tools/map.html))
7. Creates tab-separated count files.

## Installation and basic usage

You will need to install [`Conda`](https://docs.conda.io/en/latest/) package manager.

1. Open terminal and clone repository
```
git clone https://github.com/baldikacti/tnseq-homeostasis.git
```

2. Change working directory to Preprocessing
```
cd tnseq-homeostasis/Preprocessing
```

3. Open `config/config.yaml` in your favorite text editor. 

Change the `fastq: "data/test/"` paramater to the directory that contains your fastq files

Change the `results: "results/test/"` paramater to the directory that you want the results to be exported to.

It is important to keep the forward slashed at the end of directory paths.

The config file contains parameters for fasta and bed formatted gene feature references. Change those as needed for different organisms.

4. Run the pipeline with select number of cores.

```
snakemake --cores 8
```

## Outputs

1. **preprocess** directory containing filtered fastq files based on the transposon sequence
2. **bwa_aln** directory containing BAM files.
3. **bigwig** directory containing BigWig files.
4. **mapped** directory containing mapped bed formatted files
5. **read_counts** directory containing:
   1. `totalcounts_mid`: total count files that are mapped to clipped gene features.
   2. `totalcounts_full`: total count files that are mapped to full sized gene features.
   3. `uniquecounts_mid`: unique count files that are mapped to clipped gene features.
   4. `uniquecounts_full`: unique count files that are mapped to full sized gene features.
