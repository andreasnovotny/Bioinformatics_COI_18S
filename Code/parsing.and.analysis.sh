# High throuput processing of multiple MiSeq runs to be equaly treated.
# This script contains all the file parsing, and executions to run DADA2
# analysis of samples associated to the QU39 time series. 18S and COI
# sequencing projects.

# Author: Andreas Novotny
# Contact: mail@andreasnovotny.se


# Preparationn
cd Bioinformatics_COI_18S
mkdir Data

# Download Metazoogene database and convert to dada2 trainset
mkdir Data/Metazoogene
Rscript Code/MZGdb_to_DADA2.R
gzip Data/Metazoogene/*


# Define function for MiSeq analysis:
analyse_MiSeq () { # $1 = Project name, $2 = source fastqfiles, $3 = COI or 18S

    # Create analsis directory and copy files
    mkdir Data/$1 Data/$1/Fastq
    cp $2/*fastq.gz Data/$1/Fastq

    # Run dada2 pipeline
    if [ $3 = 18S ]; then
        Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/Bioinformatics_COI_18S/Data/$1"
    fi

    if [ $3 = COI ]; then
        Rscript Code/processing.COI.dada2.R "/home/andreas.novotny/Bioinformatics_COI_18S/Data/$1"
    fi

    # Copy files to provessed data
    mkdir ProcessedData/$1
    cp Data/$1/ASV/sequence_ASVname_mapping.txt ProcessedData/$1
    cp Data/$1/ASV/sequence_table.merged.w_ASV_names.txt ProcessedData/$1
    cp Data/$1/ASV/sequence_ASVname_mapping.txt ProcessedData/$1
    cp Data/$1/Taxonomy/* ProcessedData/$1

    # Remove all fastq files to clear up space
    rm -r Data/$1/Fastq
}


# Example executions - CHANGE ME:
analyse_MiSeq name_of_MiSeq_project /full/path/to/archive/containing/fastq/files Gene

# 18S MiSeq data analysis
###############################
# analyse_MiSeq 18S_QPKbulk_2017 /mnt/Genomics/MiSeq/18S_QPKbulk_2017_Run20230331/Run20230331/Alignment_1/20230403_003911/Fastq 18S
# keep adding analyses....

# COI analysis
################################
# analyse_MiSeq COI_QU39-2017 /mnt/Genomics/MiSeq/COI_Q39_2017_Run20220414/Run20220414/Alignment_1/20220417_000816/Fastq COI
# keep adding analyses....

### Expected Output:

# MiSeq Libraries are analyzed one at the time.
# 1. For each MiSeq run a called name_of_MiSeq_project will be created in ./Data
# 2. ./Data/name_of_MiSeq_project/FASTQ will contain a copy of raw sequences,
#   as well as trimmed sequences.
# 3. ./Data/name_of_MiSeq_project/Report will contain report metrics, tables and figures,
#   and will be updated as the script runs.
# 4. ./Data/name_of_MiSeq_project/ASV will contain ASV tables, taxonomy tables, and ASV name to sequence maping.
# 5. At the end of the script, the FASTQ folder will be deleated. This allows for
#   analyzing multiple runs at the same time without running out of space.
#   The content of the ASV folder will be copied to ./ProcessedData/name_of_MiSeq_project/.