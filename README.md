**Assembly and quantification metatranscriptome using metagenome data**.

Version: see VERSION

## Introduction

**MetaGT** is a bioinformatics analysis pipeline used for improving and quantification 
metatranscriptome assembly using metagenome data. The pipeline supports Illumina sequencing 
data and complete metagenome and metatranscriptome assemblies. The pipeline involves the 
alignment of metatranscriprome assembly to the metagenome assembly with further extracting CDSs,
which are covered by transcripts.

The pipeline is built using Nextflow, a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible. The Nextflow DSL2 implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install any of [`Conda`](https://conda.io/miniconda.html) for full pipeline reproducibility 

3. Download the pipeline, e.g. by cloning metaGT GitHub repository:

    ```bash
    git clone git@github.com:ablab/metaGT.git
    ```
   
4. Test it on a minimal dataset by running:

    ```bash
    nextflow run metaGT -profile test,conda
    ```
   
5. Start running your own analysis!
    > Typical command for analysis using reads:

    ```bash
    nextflow run metaGT -profile <conda> --dna_reads '*_R{1,2}.fastq.gz' --rna_reads '*_R{1,2}.fastq.gz'
    ```
    > Typical command for analysis using multiple files with reads:

    ```bash
    nextflow run metaGT -profile <conda> --dna_reads '*.yaml' --rna_reads '*.yaml' --yaml
    ```
    > Typical command for analysis using assemblies:

    ```bash
    nextflow run metaGT -profile <conda> --genome '*.fasta' --transcriptome '*.fasta'
    ```
## Pipeline Summary
Optionally, if raw reades are used:

<!-- TODO nf-core: Fill in short bullet-pointed list of default steps of pipeline -->

* Sequencing quality control (`FastQC`)
* Assembly metagenome or metatranscriptome (`metaSPAdes, rnaSPAdes `)

By default, the pipeline currently performs the following:

* Annotation metagenome (`Prokka`)
* Aligning metatranscriptome on metagenome (`minimap2`)
* Annotation unaligned transcripts (`TransDecoder`)
* Clustering covered CDS and CDS from unaligned transcripts (`MMseqs2`)
* Quantifying abundances of transcripts (`kallisto`)

## Citation

MetaGT was developed by Daria Shafranskaya and Andrey Prjibelski.
If you use it in your research please cite:

[MetaGT: A pipeline for de novo assembly of metatranscriptomes with the aid of metagenomic data](https://doi.org/10.3389/fmicb.2022.981458)

## Feedback and bug report

If you have any questions, please leave an issue at out [GitHub page](https://github.com/ablab/metaGT/issues).
