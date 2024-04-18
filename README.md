# AmpliconSeqAnalysis
scripts for our data processing pipelines

## Order of operation

1. Code/analyze.all.datasets.sh

This script is used to copy raw data from MiSeq folder an run bioinformatic analysis of 12S, COI, and 18S respectively.
The input s and outputs of all these processes will end up in the Data directory, that has one subdirectory per sequencing run.
The entire Data directory is under gitignore.
The final outputs of 1. will be coppied to Processed data.


2. prepare_datasets.Rmd

This interactive script is used to streamline and merge all datasets into a final sharable format.
The data is packed into standard formats including one ASV, Taxtab, and Samptab per gene.
The final products of this script is copied to Google Drive, to be used in downstream analyses.