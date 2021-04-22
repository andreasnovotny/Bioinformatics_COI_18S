#code for producing taxonomic assignments using blast, for CO1 and 12S amplicon data
#author: Evan Morien
#last modified: April 21st, 2021

#Introduction
#	The following code/pipeline was developed over the course of 2020/2021 to use blast to generate taxonomic assignments for both CO1 and 12S amplicon sequencing experiments
#	Code is divided into sections based on the amplicon in question, as well as the reference database

####12S amplicon####

#assign taxonomy for dada2-processed 12S amplicon data with blast using the NCBI NT database
#RDP assignments for 12S should be preferred, but this supplementary method has the added benefit of identifying bacterial and human contaminant sequences
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2020_08_28/blastdb/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query 12S_ASV_sequences.length_var.fasta  -out blast_96_sim/12S_ASV_sequences.length_var.blast.out
python2 ~/programs/Simple-LCA/add_taxonomy.w_superkingdom.py -i blast_96_sim/12S_ASV_sequences.length_var.blast.out -t ~/programs/Simple-LCA/rankedlineage.dmp -m ~/programs/Simple-LCA/merged.dmp -o blast_96_sim/12S_ASV_sequences.length_var.taxonomy_added.out
python2 ~/programs/Simple-LCA/lca.w_superkingdom.py -i blast_96_sim/12S_ASV_sequences.length_var.taxonomy_added.out -o taxonomy_table.12S.NCBI_NT.96sim.txt -b 100 -id 96 -cov 50 -t yes -tid 98 -tcov 50 -fh environmental,unidentified,phylum -flh unclassified


####CO1 amplicon####

#assign taxonomy for dada2-processed CO1 amplicon data with blast using our Hakai CO1 barcode sequences
####taxonomy assignment with blast alone, using hakai barcode blast DB####
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 1 -perc_identity 97 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/CO1_database/hakai_barcode_blast_DB/hakai_barcode_blast_DB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.fasta  -out taxonomy_table.CO1.merged.hakai_barcodeDB.97sim_cutoff.txt


#assign taxonomy for dada2-processed CO1 amplicon data with blast using NCBI NT and custom CO1 databases together
# blast against custom blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/CO1_database/blast_DB/CO1.BOLD_genbank_combined.rep_set.blast_DB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.fasta  -out CO1_ASV_sequences.customDB.blast.out
# blast against genbank NT blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2020_08_28/blastdb/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.fasta  -out CO1_ASV_sequences.blast.out

#IMPORTANT: next steps are done in R, for simplicity's sake. a custom script here would also work. this is simpler.
#we need to add taxonIDs for the customDB (adding them directly to the blast DB has not worked in the past, they don't get returned in the blast output). Using the blast output and a map of accessions to taxonIDs, we add the taxonIDs for each blast result.
library(tidyverse)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)
taxonidmap <- read.delim("~/projects/taxonomyDBs/CO1_database/taxonID_map/CO1.BOLD_genbank_combined.taxonID_map.w_hakai_barcodes.txt", sep="\t", header=F)
blastfile <- read.delim("CO1_ASV_sequences.customDB.blast.out", sep="\t", header=F)
colnames(taxonidmap) <- c("accession", "taxonID")
colnames(blastfile) <- c("asv", "col2", "accession", "blasttaxid", "col5", "col6", "col7", "col8")
taxonidmap$accession <- trimws(taxonidmap$accession, which = c("both"))
blastfile_wtaxids <- merge(blastfile,taxonidmap, by="accession", all.x=TRUE)
blastfile_output <- subset(blastfile_wtaxids, select=c("asv", "col2", "accession", "taxonID", "col5", "col6", "col7", "col8")) #it's okay here that "col2" is just all NA values, we need it to conform to the blast output from the NT database, but it doesn't get used by the taxonomy assignment scripts
blastfile_output <- blastfile_output[order(blastfile_output$asv),]
write.table(blastfile_output, "CO1_ASV_sequences.customDB.blast.out", row.names=F, col.names=F, quote=F, sep="\t") #overwriting input

#combine blast results in a way that the add taxonomy and LCA scripts can handle
blastout_customDB <- read.delim("CO1_ASV_sequences.customDB.blast.out", sep="\t", header=F)
blastout_NCBINT <- read.delim("CO1_ASV_sequences.blast.out", sep="\t", header=F)
blastout_combined <- rbind(blastout_customDB, blastout_NCBINT) #combine tables
tmp <- blastout_combined[order(-blastout_combined$V5),] #order descending order for percent identity
tmp <- tmp[order(tmp$V1),] #order by ASV name
blastout_combined <- tmp 
#write to file
write.table(blastout_combined, "CO1_ASV_sequences.combined.blast.out", row.names=F, col.names=F, quote=F, sep="\t")
#now quit R and continue with the remaining code in your bash shell

#avoid " keyerror: 'NA' " with python script by filtering on "NA" as a wholeword string
#explanation: occasionally, an output line from blast has an NA, which causes an error with the Simple-LCA script below. quick fix is to remove these lines (they're quite rare anyway)
grep -v -w "NA" CO1_ASV_sequences.combined.blast.out > tmp

#execute first step for the LCA program (adding taxonomy strings based on taxonIDs in blast output)
python2 ~/programs/Simple-LCA/add_taxonomy.w_superkingdom.py -i tmp -t ~/programs/Simple-LCA/rankedlineage.dmp -m ~/programs/Simple-LCA/merged.dmp -o CO1_ASV_sequences.taxonomy_added.out
rm tmp #remove intermediate file we created earlier

#in this section, a series of taxonomy string modifications are made. This makes the final output more readable/easier to understand, but is optional.
#IMPORTANT: please note that the filtering criteria in the final step depend on some of this filtering (e.g. blast hits with the word "phylum" will be removed. see -fh parameter in final step) 
#	because i'm replacing "unknown phylum" in these sed commands, these sequences are retained.
#	if you choose not to do the replacement, the blast hits with "unknown plylum" will not be used in LCA determination unless you also change the filtering criteria set in the -fh parameter during the final step.
#	also note, "unknown phylum" is present in any taxonomy where the clade's phylum is uncertain in the NCBI taxonomy system, it doesn't indicate any other kind of uncertainty about the provenance of the sequence.

#label fix for clades missing "kingdom" label
sed -i 's/unknown kingdom \/ Bacillariophyta/Bacillariophyta \/ Bacillariophyta/g' CO1_ASV_sequences.taxonomy_added.out #Bacillariophyta
sed -i 's/unknown kingdom \/ Ciliophora/Ciliophora \/ Ciliophora/g' CO1_ASV_sequences.taxonomy_added.out #Ciliophora
sed -i 's/unknown kingdom \/ Discosea/Discosea \/ Discosea/g' CO1_ASV_sequences.taxonomy_added.out #Discosea
sed -i 's/unknown kingdom \/ Evosea/Evosea \/ Evosea/g' CO1_ASV_sequences.taxonomy_added.out #Evosea
sed -i 's/unknown kingdom \/ Haptista/Haptista \/ Haptista/g' CO1_ASV_sequences.taxonomy_added.out #Haptista
sed -i 's/unknown kingdom \/ Rhodophyta/Rhodophyta \/ Rhodophyta/g' CO1_ASV_sequences.taxonomy_added.out #Rhodophyta

#and for those missing kingdom + phylum labels
sed -i 's/unknown kingdom \/ unknown phylum \/ Chrysophyceae/Chrysophyceae \/ Chrysophyceae \/ Chrysophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Chrysophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Cryptophyceae/Cryptophyceae \/ Cryptophyceae \/ Cryptophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Cryptophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Oomycota/Oomycota \/ Oomycota \/ Oomycota/g' CO1_ASV_sequences.taxonomy_added.out #Oomycota
sed -i 's/unknown kingdom \/ unknown phylum \/ Phaeophyceae/Phaeophyceae \/ Phaeophyceae \/ Phaeophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Phaeophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Phaeophyceae/Phaeophyceae \/ Phaeophyceae \/ Phaeophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Phaeophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Bigyra/Bigyra \/ Bigyra \/ Bigyra/g' CO1_ASV_sequences.taxonomy_added.out #Bigyra
sed -i 's/unknown kingdom \/ unknown phylum \/ Dictyochophyceae/Dictyochophyceae \/ Dictyochophyceae \/ Dictyochophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Dictyochophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Dinophyceae/Dinophyceae \/ Dinophyceae \/ Dinophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Dinophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Pelagophyceae/Pelagophyceae \/ Pelagophyceae \/ Pelagophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Pelagophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Raphidophyceae/Raphidophyceae \/ Raphidophyceae \/ Raphidophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Raphidophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Synurophyceae/Synurophyceae \/ Synurophyceae \/ Synurophyceae/g' CO1_ASV_sequences.taxonomy_added.out #Synurophyceae

#label for those missing kindom + phylum + class labels
sed -i 's/unknown kingdom \/ unknown phylum \/ unknown class \/ Telonemida/Telonemida \/ Telonemida \/ Telonemida \/ Telonemida/g' CO1_ASV_sequences.taxonomy_added.out #Telonemida
sed -i 's/unknown kingdom \/ unknown phylum \/ unknown class \/ Jakobida/Jakobida \/ Jakobida \/ Jakobida \/ Jakobida/g' CO1_ASV_sequences.taxonomy_added.out #Jakobida

#label fix for bacterial assignments
sed -i 's/Bacteria \/ unknown kingdom/Bacteria \/ Bacteria/g' CO1_ASV_sequences.taxonomy_added.out #Bacteria have no "kingdom" label in NCBI, only superkingdom and then phylum. adding "bacteria" as a kingdom label in our taxonomy assignments makes them more informative

#execute final step for the LCA program (forming consensus taxonomies for each ASV)
python2 ~/programs/Simple-LCA/lca.w_superkingdom.py -i CO1_ASV_sequences.taxonomy_added.out -o taxonomy_table.CO1.merged.combined_DB.96sim_blast_LCA.txt -b 100 -id 96 -cov 50 -t yes -tid 98 -tcov 50 -fh environmental,unidentified,phylum -flh unclassified
