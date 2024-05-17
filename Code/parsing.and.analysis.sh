# High throuput processing of multiple MiSeq runs to be equaly treated.
# This script contains all the file parsing, and executions to run DADA2
# analysis of samples associated to the QU39 time series. 18S and COI
# sequencing projects.

# Author: Andreas Novotny
# Contact: mail@andreasnovotny.se


# Preparationn
cd AmpliconSeqAnalysis
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
        Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/$1"
    fi

    if [ $3 = COI ]; then
        Rscript Code/processing.COI.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/$1"
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

# 18S MiSeq data analysis
###############################
analyse_MiSeq 18S_QPKbulk_2017 /mnt/Genomics/MiSeq/18S_QPKbulk_2017_Run20230331/Run20230331/Alignment_1/20230403_003911/Fastq 18S
analyse_MiSeq 18S_QU39_1 /mnt/RawData/QU39-18S-data/1stPool_Plates1-2-3/rerun/KelloggHakai18SPool1Rerun 18S
analyse_MiSeq 18S_QU39_2 /mnt/RawData/QU39-18S-data/2ndPool_Plates10-11-12/full-run/Hakai18SPool2_finalV3/ 18S
analyse_MiSeq 18S_QU39_3 /mnt/RawData/QU39-18S-data/3rdPool_plates4-5-6/full-run/Hakai18SPool3_finalV3 18S
analyse_MiSeq 18S_QU39_4 /mnt/RawData/QU39-18S-data/HakaiEUKpool4/ 18S
analyse_MiSeq 18S_QU39_5 /mnt/RawData/QU39-18S-data/HakaiEUKpool5 18S
analyse_MiSeq 18S_QU39_6 /mnt/Genomics/MiSeq/18S_HOME6_Run20240301/Run20240301/Alignment_1/20240304_025441/Fastq 18S

# COI analysis
################################
analyse_MiSeq COI_QPKbulk_2017 /mnt/Genomics/MiSeq/COI_QPKbul_Run20230321/COI_requeued_Mar31/230321_M06773_0116_000000000-KR4YK/Alignment_3/20230331_132717/Fastq COI
analyse_MiSeq COI_QU39-2017 /mnt/Genomics/MiSeq/COI_Q39_2017_Run20220414/Run20220414/Alignment_1/20220417_000816/Fastq COI
analyse_MiSeq COI_QU39-2018 /mnt/Genomics/MiSeq/COI_Q392018_Run20230203/Run20230203/Alignment_1/20230205_212315/Fastq COI
analyse_MiSeq COI_QU39-2019 /mnt/Genomics/MiSeq/COI_Q392019_Run20230209/Run20230209/Alignment_1/20230211_215800/Fastq COI
analyse_MiSeq COI_Zoopsprint2022 /mnt/Genomics/MiSeq/COI_Zoosprint2022_Run20221125/Run20221125/Alignment_1/20221127_181852/Fastq COI

# Merge datasets and prepare to export
################################

# Continue to prepare_datasets.Rmd (currently interractive)