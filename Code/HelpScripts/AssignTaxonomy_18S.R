

library(dada2)
library(tidyverse)
seqs <- read.delim("ASV/sequence_table.merged.txt", sep = " ") %>% 
  pivot_longer(2:length(.), names_to = "Sample", values_to = "abundance") %>% 
  rename(sequence = Sequence) %>%
  as.data.frame()

taxa_boot <- assignTaxonomy(seqs,
    "/data/taxonomyDBs/PR2/pr2_version_4.14.0_SSU_dada2.fasta.gz",
    multithread=TRUE,
    taxLevels = c("Kingdom","Supergroup","Division","Class","Order","Family","Genus","Species"),
    outputBootstraps = TRUE)

saveRDS(taxa_boot, "tax_tab_18S.RDS")