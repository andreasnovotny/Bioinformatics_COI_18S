
wd='/home/andreas.novotny/AmpliconSeqAnalysis'
cd ${wd}

# COI analysis
################################
cp -r /mnt/MiSeq/COI_Q39_2017_Run20220414/Run20220414/Alignment_1/20220417_000816/Fastq Data/COI_QU39-2017
Rscript Code/processing.CO1.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/COI_QU39-2017"

cp -r /mnt/MiSeq/COI_Q392018_Run20230203/Run20230203/Alignment_1/20230205_212315/Fastq Data/QU39-2018
Rscript Code/processing.CO1.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/COI_QU39-2018"

cp -r /mnt/MiSeq/COI_Q392019_Run20230209/Run20230209/Alignment_1/20230211_215800/Fastq
Rscript Code/processing.CO1.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/COI_QU39-2019"

cp -r /mnt/MiSeq/COI_Zoosprint2022_Run20221125/Run20221125/Alignment_1/20221127_181852/Fastq
Rscript Code/processing.CO1.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/COI_Zoopsprint2022"

cp -r /mnt/MiSeq/COI_QPKbul_Run20230321/COI_requeued_Mar31/230321_M06773_0116_000000000-KR4YK/Alignment_3/20230331_132717/Fastq/
Rscript Code/processing.CO1.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/COI_QPKbulk_2017"

sh Code/assign.CO1.taxonomy.sh

# 18S analysis
###############################
mkdir Data/18S_QU39_1/Fastq Data/18S_QU39_2/Fastq Data/18S_QU39_3/Fastq Data/18S_QU39_4/Fastq
mkdir Data/18S_QU39_5/Fastq
cp /mnt/RawData/QU39-18S-data/1stPool_Plates1-2-3/rerun/KelloggHakai18SPool1Rerun/*fastq.gz Data/18S_QU39_1/Fastq
cp /mnt/RawData/QU39-18S-data/2ndPool_Plates10-11-12/full-run/Hakai18SPool2_finalV3/*fastq.gz Data/18S_QU39_2/Fastq
cp /mnt/RawData/QU39-18S-data/3rdPool_plates4-5-6/full-run/Hakai18SPool3_finalV3/*fastq.gz Data/18S_QU39_3/Fastq
cp /mnt/RawData/QU39-18S-data/HOME-18S-run4/*fastq.gz Data/18S_QU39_4/Fastq
cp /mnt/RawData/QU39-18S-data/HOME-18S-run5/*fastq.gz Data/18S_QU39_5/Fastq

cp -r /mnt/MiSeq/18S_QPKbulk_2017_Run20230331/Run20230331/Alignment_1/20230403_003911/Fastq Data/

Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/18S_QU39_1"
Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/18S_QU39_2"
Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/18S_QU39_3"
Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/18S_QU39_4"
Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/18S_QU39_5"


# Copy Output 
mkdir ProcessedData/18S_Andreas/18S_QU39_1
cp Data/18S_QU39_1/ASV/sequence_table.merged.txt ProcessedData/18S_Andreas/18S_QU39_1
cp Data/18S_QU39_1/tax_tab_18S_pr2.RDS ProcessedData/18S_Andreas/18S_QU39_1

mkdir ProcessedData/18S_Andreas/18S_QU39_2
cp Data/18S_QU39_2/ASV/sequence_table.merged.txt ProcessedData/18S_Andreas/18S_QU39_2
cp Data/18S_QU39_2/tax_tab_18S_pr2.RDS ProcessedData/18S_Andreas/18S_QU39_2

mkdir ProcessedData/18S_Andreas/18S_QU39_3
cp Data/18S_QU39_3/ASV/sequence_table.merged.txt ProcessedData/18S_Andreas/18S_QU39_3
cp Data/18S_QU39_3/tax_tab_18S_pr2.RDS ProcessedData/18S_Andreas/18S_QU39_3

mkdir ProcessedData/18S_Andreas/18S_QU39_4
cp Data/18S_QU39_4/ASV/sequence_table.merged.txt ProcessedData/18S_Andreas/18S_QU39_4
cp Data/18S_QU39_4/tax_tab_18S_pr2.RDS ProcessedData/18S_Andreas/18S_QU39_4

mkdir ProcessedData/18S_Andreas/18S_QU39_5
cp Data/18S_QU39_5/ASV/sequence_table.merged.txt ProcessedData/18S_Andreas/18S_QU39_5
cp Data/18S_QU39_5/tax_tab_18S_pr2.RDS ProcessedData/18S_Andreas/18S_QU39_5

mkdir ProcessedData/18S_Andreas/18S_QPKbulk_2017
cp Data/18S_QPKbulk_2017/ASV/sequence_table.merged.txt ProcessedData/18S_Andreas/18S_QPKbulk_2017
cp Data/18S_QPKbulk_2017/tax_tab_18S_pr2.RDS ProcessedData/18S_Andreas/18S_QPKbulk_2017

# 12S analysis
###############################
cp -r /mnt/MiSeq/12S_Q39_2018_Run20230329/Run20230329/Alignment_1/20230331_065301/Fastq
Rscript Code/processing.12S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/12S_QU39_2018" "MiFish-E"

cp -r /mnt/MiSeq/12S_Q39_2016_Run20230804/Run209230804/Alignment_1/20230806_011559/Fastq Data/12S_QU39_2016
Rscript Code/processing.12S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/12S_QU39_2016" "MiFish-U"

cp -r /mnt/MiSeq/12S_Q39_2017_Run20230814/Run20230814/Alignment_1/20230816_055824/Fastq Data/12S_QU39_2017
Rscript Code/processing.12S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/12S_QU39_2017" "MiFish-U"

sh Code/assign.12S.taxonomy.sh