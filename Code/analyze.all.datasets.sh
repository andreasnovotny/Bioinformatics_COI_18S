
wd='/home/andreas.novotny/AmpliconSeqAnalysis'
cd ${wd}

# COI analysis
################################

## 
cp -r /mnt/MiSeq/COI_Q39_2017_Run20220414/Run20220414/Alignment_1/20220417_000816/Fastq Data/COI_QU39-2017
Rscript Code/processing.CO1.dada2.R

cp -r 


cp -r 


cp -r 


cp -r 


sh Code/assign.CO1.taxonomy.sh


# 18S analysis
cp -r
Rscript Code/processing.18S.dada2.R "/home/andreas.novotny/AmpliconSeqAnalysis/Data/18S_QPKbulk_2017"
sh Code/assign.18S.taxonomy.sh


# 12S analysis
cp -r
Rscript 