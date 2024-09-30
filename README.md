# Metabarcoding bioinformatic analysis of 18S and COI amplicon sequences using DADA2 with the MetaZooGene and PR2 databases.

This script contains code and processed data for COI and 18S sequencing projects form the QU39 time series. (see <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1141475> for background and metadata).

## Method:

This pipeline is optimized for Illumina MiSeq libraries prepared using the following two protocols:

1.  COI - Leray primers: <https://dx.doi.org/10.17504/protocols.io.rm7vzjjjxlx1/v1>
2.  18S - Balzano primers: <https://dx.doi.org/10.17504/protocols.io.261ge5d7yg47/v1>

Primer sequences are removed from the demultiplexed reads using cutadapt version 3.7 (Martin, 2011), and any sequences with missing or erroneous primers are removed. Filtering, dereplication, error rate modeling, inference of amplicon sequence variants (ASVs), chimera detection and removal, and taxonomic assignment are accomplished using the dada2 R package version 1.30 (Callahan et al., 2016). In the filterAndTrim function, forward reads are truncated after 220 nucleotides, and reverse reads are truncated after 200 nucleotides. Sequences with an error rate (EE) higher than 3 (forward reads) and 4 (reverse reads), or sequences containing ambiguous bases are removed. Paired end reads with less than 12 nucleotide overlaps are also removed. ASVs only observed once are removed before taxonomic assignment. ASVs are annotated using the MetaZooGene database v3.0 (Bucklin et al., 2021) for COI and 18S amplicons separately. A subset of the database including all marine flora, fauna, and microbes with known appearance in the North Pacific Ocean (mode-A) is used.

## Directory and order of operation:

### ./Code

This directory contains all scripts for the bioinformatic analysis:

1.  **MZGdb_to_DADA2.R** is used to download the MetaZooGene database and translate it into a DADA2 train set. This should be executed before the other scripts. Execute like this:

``` bash
mkdir Data Data/Metazoogene
Rscript Code/MZGdb_to_DADA2.R
gzip Data/Metazoogene/*
```

2.  **processing.COI.dada2.R** and **processing.18S.dada2.R** are the scripts for running the full bioinformatic analysis for COI and 18S respectively. The scripts are optimized to be executed for one MiSeq sequencing library at the time. The output of the scipt is stored in ./ProcessedData. An example execution:

``` bash
mkdir Data/Example Data/Example/Fastq
cp path_to_raw_data/*fastq.gz Data/Example/Fastq

Rscript Code/processing.COI.dada2.R "/full/path/to/Bioinformatics_COI_18S/Data/Example"
```

3.  **parsing.and.analysis.sh** is the working script to parse files and reproduce the analysis of the Hakai QU39 datasets. This script is used to copy raw data from MiSeq folder an run bioinformatic analysis of 12S, COI, and 18S respectively. The input s and outputs of all these processes will end up in the Data directory, that has one subdirectory per sequencing run. The entire Data directory is under gitignore. The final outputs of 1. will be coppied to Processed data.

### ./Code/HelpScripts

This folder contains scripts used to support file parsing in the processing scripts above.

### ./ProcessedData

Contains the output of the DADA2 pipelines, executed by #2 above.

### ./ (home)

4.  **prepare_datasets.Rmd** interactive script used to streamline and merge all datasets into a final sharable format. The data is packed into standard formats including one ASV, Taxtab, and Samptab per gene. The final products of this script is copied to Google Drive, to be used in downstream analyses.
