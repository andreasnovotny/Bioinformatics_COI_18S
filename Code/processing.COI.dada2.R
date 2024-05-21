#!/usr/bin/Rscript

#pipeline for processing COI amplicon sequencing data
#author: Evan Morien
#modified by: Andreas Novotny


###########################################################
# Execute from shell script or command line:
# Rscript processing.COI.dada2.R "/path/to/data/directory" "COI"
#
# Revisit lines marked CHANGE ME before executing script
###########################################################

#### Intro ####
# The pipleline assumes a starting point of a project directory with a
# subfolder containing raw sequencing data the raw data may be in paired
# end (standard pipline below) or non-paired (skip steps for reverse reads,
# code will need modification to exclude references to reverse reads) format,
# but each sample should have a forward and reverse (if applicable) read
# file associated


#########################################
#### Section 1: Prepare Environment #####
#########################################

#### libraries ####

# These are the libraries used in this pipeline. in a few cases,
# the order of library loading matters, so best not to modify it.

library(dada2)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)
theme_set(theme_bw())


#### File Path Setup ####
wd <- commandArgs(TRUE)[1] #should contain analysis directory.
setwd(wd)

# CHANGE ME to sequence Input
path <- file.path(wd, "Fastq")

# Create Cutadapt output
path.cut <- file.path(path, "cutadapt")
if (!dir.exists(path.cut)) dir.create(path.cut)

# Create Report directory
path.report <- file.path(wd, "Report/")
if (!dir.exists(path.report)) dir.create(path.report)

# Prepare a ASV report file:
appendASV <- function(...) {
  cat(...,
      file = "Report/Progress_report.txt",
      sep = "\t", append = TRUE)
}

appendASV(" \n Report for ASV analysis:", wd,
          "\n Date:", date()
)

# Create taxonomy result folder:
path.tax <- file.path(wd, "Taxonomy/")
if (!dir.exists(path.tax)) dir.create(path.tax)

# Create ASV result file:
path.ASV <- file.path(wd, "ASV/")
if (!dir.exists(path.ASV)) dir.create(path.ASV)

#CHANGE ME to match the pattern of all your R1 files
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# change the delimiter in quotes and the number at the end of this command to
# decide how to split up the file name, and which element to extract for a
# unique sample name
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# COI primers
FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"
REV <- "TANACYTCNGGRTGNCCRAARAAYCA"

appendASV(" \n FWD Primers:", FWD,
          "\n REV Primers:", REV)

            
# Create Cutadapt output
path.cut <- file.path(path, "cutadapt")
if (!dir.exists(path.cut)) dir.create(path.cut)


##################################################
#### Section 2: Cut Adapters, Filter and Trim ####
##################################################

#### Plot Quality Scores #####

#randomly select a set of 49 samples
a <- sample(fnFs, ifelse(length(fnFs) < 49, length(fnFs), 49))

#identify the indices of those samples
b <- which(fnFs %in% a)
# use the indices to create two lists of corresponding fwd and rev files to plot
plotfnFs <- fnFs[b]
plotfnRs <- fnRs[b]

appendASV(" \n \n - Plotting Quality scores")

# this plots the quality profiles for each sample
pdf("Report/quality_plots.dada2.R1s.pdf", width = 16, height = 9)
  plotQualityProfile(plotfnFs) #
dev.off()
pdf("Report/quality_plots.dada2.R2s.pdf", width = 16, height = 9)
  plotQualityProfile(plotfnRs)
dev.off()

# Create all orientations of the input sequence
# The Biostrings works w/ DNAString objects rather than character vectors
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna),
                Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector

} # End of allOrients function

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Put N-filterd files in filtN/ subdirectory
fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

appendASV(" \n - Removing ambigous sequences")
# This initial filter and trim will only remove sequences that has "Ns".
# This is to improve accuracy of cutadapt.
tmp_out <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN,
              trimLeft = c(0, 0),
              maxN = 0,
              multithread = 32,
              compress = TRUE,
              matchIDs = TRUE)


retained <- as.data.frame(tmp_out)
retained$percentage_retained <- retained$reads.out / retained$reads.in * 100

write.table(retained,
            "Report/retained_reads.filterAndTrim_step.csv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
} # End of primerHits function



index <- ifelse(length(fnFs) < 5, length(fnFs), 5)
# this is the index of the file we want to check for primers,
# within the lists "fn*s.filtN", it can be any number from 1
# to N, where N is the number of samples you are processing

# Construct table of detected primers.
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits,
                                fn = fnFs.filtN[[index]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits,
                                fn = fnRs.filtN[[index]]),
      REV.ForwardReads = sapply(REV.orients, primerHits,
                                fn = fnFs.filtN[[index]]),
      REV.ReverseReads = sapply(REV.orients, primerHits,
                                fn = fnRs.filtN[[index]])) %>%

      write.table("Report/detected_primers.csv", sep = "\t",
                  row.names = TRUE, col.names = TRUE, quote = FALSE)
                  
#### primer removal with Cutadapt ####

# CHANGE ME to the cutadapt path on your machine
cutadapt <- "/usr/local/bin/cutadapt"

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

appendASV(" \n - Running Cutadapt")
#Run Cutadapt
for (i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags,
                            "-n", 2,  # remove FWD and REV
                            "-j", 32, # sets no. threads
                            "-o", fnFs.cut[i], # output forward
                            "-p", fnRs.cut[i], # output reverse
                            fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# sanity check, should report zero for all orientations and read sets
index <- ifelse(length(fnFs) < 5, length(fnFs), 5)
# this is the index of the file we want to check for
# primers, within the lists "fn*s.cut", it can be any
# number from 1 to N, where N is the number of samples
# you are processing

rbind(FWD.ForwardReads = sapply(FWD.orients,
                                primerHits,
                                fn = fnFs.cut[[index]]),
      FWD.ReverseReads = sapply(FWD.orients,
                                primerHits,
                                fn = fnRs.cut[[index]]),
      REV.ForwardReads = sapply(REV.orients,
                                primerHits,
                                fn = fnFs.cut[[index]]),
      REV.ReverseReads = sapply(REV.orients,
                                primerHits,
                                fn = fnRs.cut[[index]])) %>%
      write.table("Report/detected_primers_cutadapt.csv", sep = "\t",
      row.names = TRUE, col.names = TRUE, quote = FALSE)

# Forward and reverse fastq filenames have the format:
#remember to change this so it matches ALL your file names!
cutFs <- sort(list.files(path.cut, pattern = "_R1_001", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001", full.names = TRUE))



#### filter and trim reads ####

#filter and trim command. dada2 can canonically handle lots of errors,
#I am typically permissive in the maxEE parameter set here, in order to
#retain the maximum number of reads possible. error correction steps built
#into the dada2 pipeline have no trouble handling data with this many expected
# errors. it is best, after primer removal, to not truncate with 18s data, or
# with data from any region in which the length is broadly variable. you may
# exclude organisms that have a shorter insert than the truncation length
# (definitely possible, good example is giardia). defining a minimum sequence
# length is best.

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

appendASV(" \n - Running Main FilterAndTrim")

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                    truncLen = c(220, 200),  # CHANGE ME
                    trimLeft = c(0, 0),      # CHANGE ME
                    trimRight = c(0, 0),     # CHANGE ME
                    minLen = c(150, 150),    # CHANGE ME
                    maxN = c(0, 0),          # CHANGE ME
                    maxEE = c(3, 4),         # CHANGE ME
                    # LOGGBOOK maxEE:
                    # with trunclen, 270, 250 maxEE:
                    #1,1 retained 25-30 percent of reads;
                    #2,2 retained 30-50 percent;
                    #3,4 retained around 60 percent
                    # with trunclen 280, 265:fwd:fwd
                    #3,4, retained 0.2 percent SIC!
                    truncQ = c(2, 2),       # CHANGE ME
                    rm.phix = TRUE, matchIDs = TRUE,
                    compress = TRUE, multithread = 32)

retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out / retained$reads.in * 100

write.table(retained,
            "Report/retained_reads.filterAndTrim_step.csv",
            sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

##################################################
#### Section 3: Sequence Dereplication ###########
##################################################

# File Path Setup
path.cut.filt <-file.path(wd, "Fastq/cutadapt/filtered")
filtFs <- sort(list.files(path.cut.filt, pattern = "_R1_001.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(path.cut.filt, pattern = "_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)


#### Learn and plot error rates ####

# The next three sections (learn error rates, dereplication, sample inference)
# are the core of dada2's sequence processing pipeline. read the dada2 paper
# and their online documentation (linked at top of this guide) for more
# information on how these steps work

appendASV(" \n - Learning Errors")
errF <- learnErrors(filtFs, multithread = 32)
errR <- learnErrors(filtRs, multithread = 32)

# assess this graph. it shows the error rates observed in your dataset. strange
# or unexpected shapes in the plot should be considered before moving on.
pdf("Report/error_rates.dada2.R1s.pdf", width = 10, height = 10)
  plotErrors(errF, nominalQ = TRUE)
dev.off()
pdf("Report/error_rates.dada2.R2s.pdf", width = 10, height = 10)
  plotErrors(errR, nominalQ = TRUE)
dev.off()


#### Sequence dereplication ####
appendASV(" \n - Sequence dereplication")
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that
# all your R objects have the same sample names in them
names(derepFs) <- sample.names
names(derepRs) <- sample.names


#### DADA sample inference ####
appendASV(" \n - DADA ASV inference")
dadaFs <- dada(derepFs, err = errF, multithread = 32)
dadaRs <- dada(derepRs, err = errR, multithread = 32)

dadaFs[[1]]
dadaRs[[1]]

#### Remove Low Sequence Samples ####
# A "subscript out of bounds" error at the next step (merging) may indicate that
# you aren't merging any reads in one or more samples. NB, NOT getting this
# error doesn't necessarily mean that all of your samples end up with more than
# 0 merged reads, as i found out while processing a large 18s dataset. your
# guess is as good as mine as to why this error does or does not appear, but
# filtering out the samples that cause it is necessary for completion of the
# pipeline.
# method used above after the filter and trim step. if you already did this but
# still got an error when merging, try the steps below:

#  track read retention, number of unique sequences after ASV inference
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(derepFs, getN),
              sapply(derepRs, getN),
              sapply(dadaFs, getN),
              sapply(dadaRs, getN))

samples_to_keep <- track[, 4] > 20 #your threshold.
# try different ones to get the lowest one that will work.
# this method accounts for dereplication/ASVs left after inference

#record names of samples you have the option of removing #be sure to note down
# what you removed for future reference
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)]

#### merge paired reads ####
#OPTION 1: version of command with no samples left out
# mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#a "subscript out of bounds" error here may indicate that you aren't merging
# any reads in one or more samples. you can remove samples with low counts
# from the workflow before the filterAndTrim step (a few steps back), or you
# can filter samples using the information from the dereplication and
# sample-inference steps (section just above)

#OPTION 2: modify command when removing low-sequence samples
appendASV(" \n - Merge Pair End reads")
mergers <- mergePairs(dadaFs[samples_to_keep],
                      derepFs[samples_to_keep],
                      dadaRs[samples_to_keep],
                      derepRs[samples_to_keep],
                      verbose = TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#### construct sequence table ####
seqtab <- makeSequenceTable(mergers)



#### View Sequence Length Distribution Post-Merging ####
# most useful with merged data. this plot will not show you much for forward
# reads only, which should have a uniform length distribution.

#tabulate sequence length distribution
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab))))
pdf("Report/length_histogram.merged_reads.pdf", width = 10, height = 8)
plot(x = length.histogram[, 1], y = length.histogram[, 2])
dev.off()


##################################################
#### Section 4: Remove Singeltons and Bimeras ####
##################################################

#### remove low-count singleton ASVs ####
appendASV(" \n \n - Remove singletons:")
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

# generate counts of sample per ASV
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV

# sanity check
# (this should be the same as last command, but the dimensions reversed)

appendASV("\n Dimentions of sequence table with singeltons:", dim(seqtab))


#create relative abundance table
otus_rel_ab <- transform_sample_counts(otus, function(x) x / sum(x))
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame

# if there are samples with no merged reads in them, and they passed the
# merge step (a possiblity, converting to a relative abundance table
# produes all NaNs for that sample. these need to be set to zero so we
# can do the calculations in the next steps.)
df[is.na(df)] <- 0

# compute row sums (sum of relative abundances per ASV. for those only
# present in one sample, this is a value we can use to filter them for
# relative abundance on a per-sample basis)
otus_rel_ab.rowsums <- rowSums(df)

# which ASVs are only present in one sample
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1)

# here is where you set your relative abundance threshold #which ASVs
# pass our filter for relative abundance
b <- which(otus_rel_ab.rowsums <= 0.001)

 #A also in B (we remove singleton ASVs that have a lower relative
 # abundance value than our threshold)
rows_to_remove <- intersect(a, b)
#filter OTU table we created earlier
otus_filt <- otus[-rows_to_remove, ]

#convert filtered OTU table back to a sequence table matrix
# to continue with dada2 pipeline
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt)))

appendASV(
  "\n dimensions of ASV table with singletons:", dim(otus),
  #how many of our singleton ASVs fail on this filter
  "\n # ASVs singeltons removed:", length(intersect(a, b)),
  #how many ASVs did you retain?
  "\n dimensions of ASV table without singletons:", dim(otus_filt))


#### remove chimeras ####
# here we remove "bimeras" or chimeras with two sources.
# look at "method" to decide which type of pooling you'd
# like to use when judging each sequence as chimeric or non-chimeric

#this step can take a few minutes to a few hours,
# depending on the size of your dataset
appendASV("\n \n - Remove Bimeras")
seqtab.nosingletons.nochim <- removeBimeraDenovo(
  seqtab.nosingletons,
  method = "pooled",
  multithread = 36,
  verbose = TRUE)

appendASV("\n dimensions of ASV table after chimera removal",
          dim(seqtab.nosingletons.nochim))

saveRDS(seqtab.nosingletons.nochim, "seqtab.nosingletons.nochim.RDS")

# proportion of nonchimeras #it should be relatively high after
# filtering out your singletons/low-count ASVs, even if you lose
# a lot of ASVs, the number of reads lost should be quite low

appendASV("\n proprtion of chimeric to non chimeric reads:",
          sum(seqtab.nosingletons.nochim) / sum(seqtab.nosingletons))

appendASV("\n attempt taxonomy assignment....")

taxa_boot <- assignTaxonomy(seqtab.nosingletons.nochim,
  "../Metazoogene/MZGdb_COI_NPac_ALL_mode-A_v3.0.fasta.gz",
  multithread=TRUE,
  taxLevels = c("Kingdom","Phylum","Subphylum","Superclass","Subsuperclass","Class",
  "Infraclass", "Superorder", "Order", "Family", "Genus", "Species"),
  outputBootstraps = TRUE)


saveRDS(taxa_boot, "Taxonomy/tax_tab_COI_MZGdb.RDS")

##################################################
#### Section 5: Finalizing quality steps #########
##################################################

appendASV("\n \n - Creating final outputs")

#### track read retention through steps ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep, ],
              sapply(dadaFs[samples_to_keep], getN),
              sapply(dadaRs[samples_to_keep], getN),
              sapply(mergers, getN),
              rowSums(seqtab.nosingletons),
              rowSums(seqtab.nosingletons.nochim))

# If processing only a single sample, remove the sapply calls: e.g.
# replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track,
              100 - track[, 6] / track[, 5] * 100,
              100 - track[, 7] / track[, 6] * 100,
              track[, 7] / track[, 1] * 100)

colnames(track) <- c(
  "input", "filtered", "denoisedF", "denoisedR", "merged",
  "nosingletons", "nochimeras", "percent_singletons",
  "percent_chimeras", "percent_retained_of_total")



#### save output from sequnce table construction steps ####
write.table(data.frame("row_names" = rownames(track), track),
                      "Report/read_retention.merged.csv",
                      row.names = FALSE, quote = FALSE, sep = "\t")

write.table(data.frame("row_names" = rownames(seqtab.nosingletons.nochim),
                      seqtab.nosingletons.nochim),
                      "ASV/sequence_table.merged.txt",
                      row.names = FALSE, quote = FALSE, sep = "\t")

#### save sequences for both ASV tables , and do taxonomy assignment ####

# save the new names and sequences as a .fasta file in your project working
# directory, and save a table that shows the mapping of sequences to new
# ASV names

# transposed (OTUs are rows) data frame. unclassing the otu_table() output
# avoids type/class errors later on
my_otu_table <- t(as.data.frame(seqtab.nosingletons.nochim))

#store sequences in character vector
ASV.seq <- as.character(unclass(row.names(my_otu_table)))
ASV.num <- paste0("ASV", seq(ASV.seq), sep = "") #create new names
write.table(cbind(ASV.num, ASV.seq),
            "ASV/sequence_ASVname_mapping.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#save sequences with new names in fasta format
write.fasta(sequences = as.list(ASV.seq),
            names = ASV.num, "ASV/ASV_sequences.fasta")

#IMPORTANT: sanity checks
testCols <- colnames(seqtab.nosingletons.nochim) == ASV.seq

#only proceed if this tests as true for all elements
appendASV("\n",
          ifelse(FALSE %in% testCols == TRUE,
                print("WARNING - ASVs and Colnames do not agree"),
                print("ASV names agree - all OK")))

#assign new ASV names
colnames(seqtab.nosingletons.nochim) <- ASV.num


#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names" = rownames(seqtab.nosingletons.nochim),
                      seqtab.nosingletons.nochim),
                      "ASV/sequence_table.merged.w_ASV_names.txt",
                      row.names = FALSE, quote = FALSE, sep = "\t")

appendASV("\n \n ASV ANALYSIS COMPLETE")