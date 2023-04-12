Rscript Code/HelpScripts/Merge_ASV.r \
    Data/Assign_18S_Taxonomy/ASV \
    Data/18S_QPKbulk_2017

wd='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_18S_Taxonomy'

cd ${wd}

library(dada2)
seqtab <- read.delim("ASV/sequence_ASVname_mapping.CO1.txt")

taxa_boot <- assignTaxonomy(seqtab[2],
    "/data/taxonomyDBs/PR2/pr2_version_4.14.0_SSU_dada2.fasta.gz",
    multithread=TRUE,
    taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
    outputBootstraps = TRUE)