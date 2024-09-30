# Code for ASV sequence inference and Taxonomic assignment using the MetaZooGene database:

Instructions for running the processing scripts. See code annotations for details.

### processing.CO1.dada2.R and processing.18S.dada2.R

Example execution:

``` bash
Rscript Code/processing.COI.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/COI_QU39-2017/"
```

The pipleline assumes a starting point of a project directory with a subfolder containing raw sequencing data the raw data may be in paired end (standard pipline below) or non-paired (skip steps for reverse reads, code will need modification to exclude references to reverse reads) format, but each sample should have a forward and reverse (if applicable) read file associated.

The scripts are divided into several subsections:

1.  Prepare environment, and parse files.

2.  Cut adapters, filter and trim sequences

    -   Plot Quality Scores
    -   Remove ambiguous sequences with filterAndTrim
    -   Primer removal with cutadapt
    -   Truncate reads, filter by FASTQ error score, once again using FilterAndTrim

3.  Sequence dereplication

    -   Learn and plot error rates
    -   Sequence dereplication with dada

4.  Bimera removal and Taxonomic assignment

5.  Track read retention and quality checks
