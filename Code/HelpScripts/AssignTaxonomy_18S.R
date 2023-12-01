

library(dada2)
library(tidyverse)
seqs <- read.delim("ASV/sequence_table.merged.txt", sep = " ") %>% 
  pivot_longer(2:length(.), names_to = "Sample", values_to = "abundance") %>% 
  rename(sequence = Sequence) %>%
  as.data.frame()


taxa_boot <- assignTaxonomy(seqs,
    "/data/taxonomyDBs/silva_for_dada2/v132_for_parfreylab/18s/silva_132.18s.99_rep_set.dada2.fa.gz",
    multithread=TRUE,
    taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species","Accession"),
    outputBootstraps = TRUE)

saveRDS(taxa_boot, "tax_tab_18S_SILVA.RDS")


taxa_boot <- assignTaxonomy(seqs,
"/data/taxonomyDBs/PR2/pr2_version_4.14.0_SSU_dada2.fasta.gz",
multithread=TRUE,
taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
outputBootstraps = TRUE)

saveRDS(taxa_boot, "tax_tab_18S_pr2.RDS")