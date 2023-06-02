
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
cp -r /mnt/MiSeq/18S_QPKbulk_2017_Run20230331/Run20230331/Alignment_1/20230403_003911/Fastq Data/
Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/18S_QPKbulk_2017"
sh Code/assign.18S.taxonomy.sh

# 12S analysis
###############################
cp -r /mnt/MiSeq/12S_Q39_2018_Run20230329/Run20230329/Alignment_1/20230331_065301/Fastq
Rscript Code/processing.12S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/12S_QU39_2018"
sh Code/assign.12S.taxonomy.sh