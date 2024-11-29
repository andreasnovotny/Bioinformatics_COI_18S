# DNA Metabarcoding Pipeline for Marine Pelagic Biodiversity Assessment - Bioinformatics

***Hakai Institute and the University of British Columbia\
\
**For Minimum Information about an Omics Protocol (MIOP) - See Below*

## Rationale

This repository contains scripts for bioinformatic analysis of fastq sequence files obtained from pair-end MiSeq HTS optimized for the following genetic markers:

1.  **COI** **(Mitochondrial Cytochrome Oxidase I) gene**\
    Targeting **Marine Invertebrates** with the primers mICOIintF / dgHCO2198R (Leray et al 2013).\
    Prepared and sequenced according to <https://dx.doi.org/10.17504/protocols.io.rm7vzjjjxlx1/v1>\
2.  **18S SSU rRNA gene**\
    Targeting **Marine Eukaryotes** with the primers Primers: V4-565F / V4981R (Balzano et al 2015) Prepared and sequenced according to <https://dx.doi.org/10.17504/protocols.io.261ge5d7yg47/v1>

Briefly, for taxonomic identification using \~320 bp section of the cytochrome oxidase I gene (COI) is amplified using the primers mICOIintF (GGWACWGGWTGAACWGTWTAYCCYCC) and dgHCO2198R (TAAACTTCAGGGTGACCAAARAAYCA), frequently used to study marine invertebrate biodiversity (Leray *et al.*, 2013). Additionally, a \~400 bp section of the 18S SSU rRNA gene (V4 region) (18S) is amplified using the primers V4-565F (CCAGCASCYGCGGTAATTCC) and V4981R (CCAGCASCYGCGGTAATTCC), designed to capture a broad diversity of marine organisms including phytoplankton, heterotrophic protists, and zooplankton (Balzano *et al.*, 2015).\
\
This pipeline assumes that samples have been multiplexed with index sequences, and sequenced on Illumina MiSeq V3 with a 2\*300 pair end setup. Raw data should have been converted to fastq format and samples de-multiplexed according to manufacturers standards. See the links above for suggested details on library preparation and sequencing parameters.

## Method Summary

These scripts follow this general pipeline for processing:

1.  **Cutadapt:** Primer sequences are removed from the demultiplexed reads using *cutadapt* version 3.7 (Martin, 2011). Here and any sequences with missing or erroneous primers were removed.

2.  **Dada2:** Filtering, dereplication, error rate modeling, inference of amplicon sequence variants (ASVs), chimera detection and removal, and taxonomic assignment were accomplished using the dada2 R package version 1.30 (Callahan *et al.*, 2016).

    2.1 In the *filterAndTrim* function, forward reads are truncated after 220 nucleotides, and reverse reads were truncated after 200 nucleotides. Sequences with an error rate (EE) higher than 3 (forward reads) and 4 (reverse reads), or sequences containing ambiguous bases are removed. Paired end reads with less than 12 nucleotide overlaps are also removed. ASVs only observed once are removed before taxonomic assignment.

3.  **Taxonomic Assignment:** ASVs are annotated using the RDP Naïve Bayesian Classifier algoritm implemented in the *AssignTaxonomy* function in the dada2 r package (Wang *et al.*, 2007; Callahan *et al.*, 2016). The training data sets for 18S and COI was constructed separately based on the **MetaZooGene database v3.0** (Bucklin *et al.*, 2021) a curated database of marine zooplankton, fish and protists, along with their observational and sampling metadata. Here, the subset provided for North Pacific Ocean species (mode-A), where only species with previously recorded presence in the North Pacific Ocean was included is downloaded and transformed (for details, see Bucklin *et al.*, 2021).

## GUIDE TO ARCHIVED METHODOLOGY {#guide-to-archived-methodology}

### Content and Order of Operation:

###  1. ./Code

This directory contains all scripts for the bioinformatic analysis:

#### **1.1. parsing.and.analysis.sh \<- *Start here***

This is the working script to parse files. This script is used to copy raw data from an archive directory and run bioinformatic analysis of COI, and 18S respectively. The inputs and outputs of all these processes will end up in the Data directory, that has one sub directory per sequencing run. The entire Data directory is under gitignore. The final outputs of 1. will be copied to Processed data.

**1.2 MZGdb_to_DADA2.R**

This script is used to download the MetaZooGene database and translate it into a DADA2 train set. This should be executed before the other scripts. Execute like this:

``` bash
mkdir Data Data/Metazoogene
Rscript Code/MZGdb_to_DADA2.R
gzip Data/Metazoogene/*
```

**1.3 processing.COI.dada2.R** and **processing.18S.dada2.R**

These are the scripts for running the full bioinformatic analysis for COI and 18S respectively. The scripts are optimized to be executed for one MiSeq sequencing library at the time. The output of the scipt is stored in ./ProcessedData. An example execution:

``` bash
mkdir Data/Example Data/Example/Fastq
cp path_to_raw_data/*fastq.gz Data/Example/Fastq

Rscript Code/processing.COI.dada2.R "/full/path/to/Bioinformatics_COI_18S/Data/Example"
```

### 2 ./ProcessedData

Will contain the outputs of the pipelines.

## Requirements

Ensure software is installed prior to executing the analysis.

Cutadapt and R.

R packages: dada2, seqinr, phyloseq, tidyverse, reshape2, stringr, data.table, broom, ape, qualpalr, ShortRead, and Biostrings.\

## Minimum Information about an Omics Protocol (MIOP)

See [MIOP_definition.md](https://github.com/BeBOP-OBON/0_protocol_collection_template/blob/main/MIOP_definition.md) for list and definitions.

+-----------------------------------+-----------------------------------------------------------------------------+
| MIOP Term                         | Value                                                                       |
+===================================+=============================================================================+
| methodology category              | Bioinformatics                                                              |
+-----------------------------------+-----------------------------------------------------------------------------+
| project                           | Biomolecular Surveys of Marine Biodiversity in the Coastal North Pacific    |
+-----------------------------------+-----------------------------------------------------------------------------+
| purpose                           | DNA metabarcoding - taxonomic assignment                                    |
+-----------------------------------+-----------------------------------------------------------------------------+
| analyses                          | Analysis of sequenced Amplicons, 18S, COI                                   |
+-----------------------------------+-----------------------------------------------------------------------------+
| geographic location               | Coastal North Pacific                                                       |
+-----------------------------------+-----------------------------------------------------------------------------+
| broad-scale environmental context | marine biome ENVO:00000447                                                  |
+-----------------------------------+-----------------------------------------------------------------------------+
| local environmental context       | coastal sea water ENVO:00002149                                             |
+-----------------------------------+-----------------------------------------------------------------------------+
| environmental medium              | sea water ENVO:00002149                                                     |
+-----------------------------------+-----------------------------------------------------------------------------+
| target                            | Marine Eukaryote Biodiversity (18S), Marine Invertebrate Biodiversity (COI) |
+-----------------------------------+-----------------------------------------------------------------------------+
| creator                           | Andreas Novotny\                                                            |
|                                   | Evan Morien                                                                 |
+-----------------------------------+-----------------------------------------------------------------------------+
| materials required                | Unix Based OS, approx. 20 parrallell CPU.                                   |
+-----------------------------------+-----------------------------------------------------------------------------+
| skills required                   | Bioinformatics                                                              |
+-----------------------------------+-----------------------------------------------------------------------------+
| time required                     | Approx. 4 h per MiSeq Library                                               |
+-----------------------------------+-----------------------------------------------------------------------------+
| personnel required                | 1                                                                           |
+-----------------------------------+-----------------------------------------------------------------------------+
| language                          | Eng. R. Python. Bash.                                                       |
+-----------------------------------+-----------------------------------------------------------------------------+
| issued                            | 2024-11-28                                                                  |
+-----------------------------------+-----------------------------------------------------------------------------+
| audience                          | scientists                                                                  |
+-----------------------------------+-----------------------------------------------------------------------------+
| publisher                         | Hakai Institute, Ocean Observing Program,\                                  |
|                                   | University of British Columbia, Pelagic Ecosystems Laboratory.              |
+-----------------------------------+-----------------------------------------------------------------------------+
| hasVersion                        | v1.0                                                                        |
+-----------------------------------+-----------------------------------------------------------------------------+
| license                           | CC BY 4.0                                                                   |
+-----------------------------------+-----------------------------------------------------------------------------+
| maturity level                    | Mature                                                                      |
+-----------------------------------+-----------------------------------------------------------------------------+

## AUTHORS

+------------------------------------------------------------------------------------------------------------------------------------+--------------------------------+
| PREPARED BY All authors known to have contributed to the preparation of this protocol, including those who filled in the template. | AFFILIATION                    |
+====================================================================================================================================+================================+
| Andreas Novotny                                                                                                                    | University of British Columbia |
+------------------------------------------------------------------------------------------------------------------------------------+--------------------------------+
| Evan Morien                                                                                                                        | Hakai Institute                |
+------------------------------------------------------------------------------------------------------------------------------------+--------------------------------+
