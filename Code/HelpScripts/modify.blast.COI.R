
# we need to add taxonIDs for the customDB
# (adding them directly to the blast DB has not worked in the past,
# they don't get returned in the blast output). Using the blast output
 # and a map of accessions to taxonIDs, we add the taxonIDs for each
 # blast result.

library(tidyverse)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)

taxonidmap <- read.delim(
    "/data/taxonomyDBs/CO1_database/taxonID_map/CO1.BOLD_genbank_combined.taxonID_map.w_hakai_barcodes.txt",
    sep = "\t", header = FALSE)

#read in blast results for 96% similarity iteration:
blastfile96 <- read.delim(
    "blast_96_sim/CO1_ASV_sequences.customDB.blast.out",
    sep = "\t", header = FALSE)

#read in blast results for 90% similarity iteration
blastfile90 <- read.delim(
    "blast_90_sim/CO1_ASV_sequences.customDB.blast.out",
    sep = "\t", header = FALSE)

#read in blast results for 90% similarity iteration
blastfile80 <- read.delim(
    "blast_80_sim/CO1_ASV_sequences.customDB.blast.out",
    sep = "\t", header = FALSE)

#join iterations for customDB blast
blastfile <- rbind(blastfile96, blastfile90, blastfile80)

colnames(taxonidmap) <-
    c("accession", "taxonID")

colnames(blastfile) <-
    c("asv", "col2", "accession", "blasttaxid", "col5", "col6", "col7", "col8")

taxonidmap$accession <- trimws(taxonidmap$accession, which = c("both"))

blastfile_wtaxids <- merge(
    blastfile, taxonidmap, by = "accession", all.x = TRUE)
#it's okay here that "col2" is just all NA values,
# we need it to conform to the blast output from the NT database,
# but it doesn't get used by the taxonomy assignment scripts
blastfile_output <- subset(
    blastfile_wtaxids,
    select = c(
        "asv", "col2", "accession", "taxonID",
        "col5", "col6", "col7", "col8"))

blastfile_output <- blastfile_output[order(blastfile_output$asv), ]

#overwriting input
write.table(
    blastfile_output,
    "blast_96_sim/CO1_ASV_sequences.customDB_96_90_80.blast.out",
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

#combine blast results in a way that the add taxonomy and LCA scripts can handle
blastout_customDB <- read.delim(
    "blast_96_sim/CO1_ASV_sequences.customDB_96_90_80.blast.out",
    sep = "\t", header = FALSE)

blastout_NCBINT96 <- read.delim(
    "blast_96_sim/CO1_ASV_sequences.blast.out",
    sep = "\t", header = FALSE)

blastout_NCBINT90 <- read.delim(
    "blast_90_sim/CO1_ASV_sequences.blast.out",
    sep = "\t", header = FALSE)

blastout_NCBINT80 <- read.delim(
    "blast_80_sim/CO1_ASV_sequences.blast.out",
    sep = "\t", header = FALSE)

#join iterations for NT blast
blastout_NCBINT <- rbind(blastout_NCBINT96, blastout_NCBINT90, blastout_NCBINT80)

 #combine tables
blastout_combined <- rbind(blastout_customDB, blastout_NCBINT)

#order descending order for percent identity
tmp <- blastout_combined[order(-blastout_combined$V5), ]

#order by ASV name
tmp <- tmp[order(tmp$V1), ]

blastout_combined <- tmp

#write to file
write.table(
    blastout_combined,
    "blast_96_sim/CO1_ASV_sequences.combined_all.blast.out",
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
#now quit R and continue with the remaining code in your bash shell