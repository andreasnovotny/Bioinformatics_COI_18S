
# Usage: Rscript Path/To/Script path/to/output path/to/inputs/multiple ...
# Use to merge ASV output from several MiSeq Runs


library(tidyverse)
library(seqinr)

readASV <- function(dirpath){

    require(tidyverse)

    dirpath %>%
        file.path("ASV/sequence_table.merged.txt") %>%
        read.delim() %>%
        column_to_rownames("row_names") %>%
        t()%>%
        as.data.frame() %>%
        rownames_to_column("Sequence")

} # End of readASV function

out_dir <- commandArgs(TRUE)[1]
in_dirs <- commandArgs(TRUE)[-1]
#in_dirs <- c("Data/DebZoop", "Data/DebZoop")
if (!dir.exists(out_dir)) dir.create(out_dir)

# Read in all ASV tables:
table_list <- mapply(readASV, in_dirs, SIMPLIFY = FALSE)

# Merge all ASV tables to one big dataframe
mergedASV <- data.frame(Sequence = as.character())

for (i in seq(table_list)){
mergedASV <- mergedASV %>%
    full_join(table_list[[i]], by = "Sequence")
}


# Save merged ASV table
write.table(mergedASV, file.path(out_dir, "sequence_table.merged.txt"))

# Rename ASVs to short names
ASV.seq <- mergedASV %>% pull(Sequence)
ASV.num <- paste0("ASV", seq(ASV.seq), sep = "") #create new names

# Save ASV name mapping
write.table(cbind(ASV.num, ASV.seq),
            file.path(out_dir, "sequence_ASVname_mapping.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#Save ASV table with ASV names
write.fasta(sequences = as.list(ASV.seq),
            names = ASV.num,
            file.path(out_dir, "ASV_sequences.fasta"))

print("Merged Sequences Rscript finalised")
