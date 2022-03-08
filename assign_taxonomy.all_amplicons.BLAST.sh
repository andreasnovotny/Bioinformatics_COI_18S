#code for producing taxonomic assignments using blast, for CO1 and 12S amplicon data
#author: Evan Morien
#last modified: June 29th, 2021

#Introduction
#	The following code/pipeline was developed over the course of 2020/2021 to use blast to generate taxonomic assignments for both CO1 and 12S amplicon sequencing experiments
#	Code is divided into sections based on the amplicon in question, as well as the reference database

#### 12S amplicon RDP with custom MitoFish classifier ####
#assign taxonomy for dada2-processed 12S amplicon data with blast using the NCBI NT database
#RDP assignments for 12S should be preferred, but this supplementary method has the added benefit of identifying bacterial and human contaminant sequences
#run 12S classifier from terrimporter on github (sequences for 12S isolated from mitofish mitochondiral genome repo, classifier trained on these sequences)
java -Xmx248g -jar ~/programs/rdp_classifier_2.13/dist/classifier.jar classify -c 0.8 -t ~/projects/taxonomyDBs/12S_database/terrimporter_12S_fish_classifier/mydata_trained/rRNAClassifier.properties -o taxonomy_table.12S.merged.RDP.txt 12S_ASV_sequences.length_var.fasta


#### 12S, 16S, or 18S blasting against NCBI NT ####
#remember to change the input and output file names so they match the amplicon you're working with
#assign taxonomy with blast NT database at 96% similarity threshold
mkdir blast_96_sim
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2020_08_28/blastdb/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query 12S_ASV_sequences.length_var.fasta  -out blast_96_sim/12S_ASV_sequences.length_var.blast.out
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim/12S_ASV_sequences.length_var.blast.out -t ~/programs/Simple-LCA/rankedlineage.dmp -m ~/programs/Simple-LCA/merged.dmp -o blast_96_sim/taxonomy
cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_12S_ASV_sequences.length_var.blast.out > tmp
python2 ~/programs/galaxy-tool-lca/lca.py -i tmp -o blast_96_sim/taxonomy_table.12S.NCBI_NT.96sim.txt -b 100 -id 96 -cov 50 -t best_hit -tid 98 -tcov 80 -fh environmental,unidentified,kingdom,phylum,class,order,family,genus -flh unclassified,unknown
#for 12S, i have found that the filtering parameters below are best for identifying bacterial (non 12S sequence)
python2 ~/programs/galaxy-tool-lca/lca.py -i tmp -o blast_96_sim/taxonomy_table.12S.NCBI_NT.96sim.txt -b 100 -id 96 -cov 50 -t best_hit -tid 98 -tcov 80 -fh environmental,unidentified,kingdom -flh unclassified


#cleanup
rm blast_96_sim/12S_ASV_sequences.length_var.blast.out #remove blast output without taxonomy
rm taxonomy_12S_ASV_sequences.length_var.blast.out #remove redundant file
mv tmp blast_96_sim/12S_ASV_sequences.length_var.blast.out #replace with taxonomy added blast output


#### CO1 amplicon blast with NT and custom CO1 sequence database ####
#assign taxonomy for dada2-processed CO1 amplicon data with blast using our Hakai CO1 barcode sequences
####taxonomy assignment with blast alone, using hakai barcode blast DB####
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 1 -perc_identity 97 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/CO1_database/hakai_barcode_blast_DB/hakai_barcode_blast_DB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.fasta  -out taxonomy_table.CO1.merged.hakai_barcodeDB.97sim_cutoff.txt

#assign taxonomy for dada2-processed CO1 amplicon data with blast using NCBI NT and custom CO1 databases together
mkdir blast_96_sim
# blast against custom blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/CO1_database/blast_DB/CO1.BOLD_genbank_combined.rep_set.blast_DB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.fasta  -out blast_96_sim/CO1_ASV_sequences.customDB.blast.out
# blast against genbank NT blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2020_08_28/blastdb/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.fasta  -out blast_96_sim/CO1_ASV_sequences.blast.out

#IMPORTANT: next steps are done in R, for simplicity's sake. a custom script here would also work. this is simpler.
#we need to add taxonIDs for the customDB (adding them directly to the blast DB has not worked in the past, they don't get returned in the blast output). Using the blast output and a map of accessions to taxonIDs, we add the taxonIDs for each blast result.
library(tidyverse)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)
taxonidmap <- read.delim("~/projects/taxonomyDBs/CO1_database/taxonID_map/CO1.BOLD_genbank_combined.taxonID_map.w_hakai_barcodes.txt", sep="\t", header=F)
blastfile <- read.delim("blast_96_sim/CO1_ASV_sequences.customDB.blast.out", sep="\t", header=F)
colnames(taxonidmap) <- c("accession", "taxonID")
colnames(blastfile) <- c("asv", "col2", "accession", "blasttaxid", "col5", "col6", "col7", "col8")
taxonidmap$accession <- trimws(taxonidmap$accession, which = c("both"))
blastfile_wtaxids <- merge(blastfile,taxonidmap, by="accession", all.x=TRUE)
blastfile_output <- subset(blastfile_wtaxids, select=c("asv", "col2", "accession", "taxonID", "col5", "col6", "col7", "col8")) #it's okay here that "col2" is just all NA values, we need it to conform to the blast output from the NT database, but it doesn't get used by the taxonomy assignment scripts
blastfile_output <- blastfile_output[order(blastfile_output$asv),]
write.table(blastfile_output, "blast_96_sim/CO1_ASV_sequences.customDB.blast.out", row.names=F, col.names=F, quote=F, sep="\t") #overwriting input

#combine blast results in a way that the add taxonomy and LCA scripts can handle
blastout_customDB <- read.delim("blast_96_sim/CO1_ASV_sequences.customDB.blast.out", sep="\t", header=F)
blastout_NCBINT <- read.delim("blast_96_sim/CO1_ASV_sequences.blast.out", sep="\t", header=F)
blastout_combined <- rbind(blastout_customDB, blastout_NCBINT) #combine tables
tmp <- blastout_combined[order(-blastout_combined$V5),] #order descending order for percent identity
tmp <- tmp[order(tmp$V1),] #order by ASV name
blastout_combined <- tmp 
#write to file
write.table(blastout_combined, "blast_96_sim/CO1_ASV_sequences.combined.blast.out", row.names=F, col.names=F, quote=F, sep="\t")
#now quit R and continue with the remaining code in your bash shell

#avoid " keyerror: 'NA' " with python script by filtering on "NA" as a wholeword string
#explanation: occasionally, an output line from blast has an NA, which causes an error with the Simple-LCA script below. quick fix is to remove these lines (they're quite rare anyway)
grep -v -w "NA" blast_96_sim/CO1_ASV_sequences.combined.blast.out > blast_96_sim/tmp

#execute first step for the LCA program (adding taxonomy strings based on taxonIDs in blast output)
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim/tmp -t ~/programs/Simple-LCA/rankedlineage.dmp -m ~/programs/Simple-LCA/merged.dmp -o taxonomy
cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_tmp > tmp

#in this section, a series of taxonomy string modifications are made. This makes the final output more readable/easier to understand, but is optional.
#IMPORTANT: please note that the filtering criteria in the final step depend on some of this filtering (e.g. blast hits with the word "phylum" will be removed. see -fh parameter in final step) 
# because i'm replacing "unknown phylum" in these sed commands, these sequences are retained.
# if you choose not to do the replacement, the blast hits with "unknown plylum" will not be used in LCA determination unless you also change the filtering criteria set in the -fh parameter during the final step.
# also note, "unknown phylum" is present in any taxonomy where the clade's phylum is uncertain in the NCBI taxonomy system, it doesn't indicate any other kind of uncertainty about the provenance of the sequence.

#label fix for clades missing "kingdom" label
sed -i 's/unknown kingdom \/ Bacillariophyta/Bacillariophyta \/ Bacillariophyta/g' tmp #Bacillariophyta
sed -i 's/unknown kingdom \/ Ciliophora/Ciliophora \/ Ciliophora/g' tmp #Ciliophora
sed -i 's/unknown kingdom \/ Discosea/Discosea \/ Discosea/g' tmp #Discosea
sed -i 's/unknown kingdom \/ Evosea/Evosea \/ Evosea/g' tmp #Evosea
sed -i 's/unknown kingdom \/ Haptista/Haptista \/ Haptista/g' tmp #Haptista
sed -i 's/unknown kingdom \/ Rhodophyta/Rhodophyta \/ Rhodophyta/g' tmp #Rhodophyta

#and for those missing kingdom + phylum labels
sed -i 's/unknown kingdom \/ unknown phylum \/ Chrysophyceae/Chrysophyceae \/ Chrysophyceae \/ Chrysophyceae/g' tmp #Chrysophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Cryptophyceae/Cryptophyceae \/ Cryptophyceae \/ Cryptophyceae/g' tmp #Cryptophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Oomycota/Oomycota \/ Oomycota \/ Oomycota/g' tmp #Oomycota
sed -i 's/unknown kingdom \/ unknown phylum \/ Phaeophyceae/Phaeophyceae \/ Phaeophyceae \/ Phaeophyceae/g' tmp #Phaeophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Phaeophyceae/Phaeophyceae \/ Phaeophyceae \/ Phaeophyceae/g' tmp #Phaeophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Bigyra/Bigyra \/ Bigyra \/ Bigyra/g' tmp #Bigyra
sed -i 's/unknown kingdom \/ unknown phylum \/ Dictyochophyceae/Dictyochophyceae \/ Dictyochophyceae \/ Dictyochophyceae/g' tmp #Dictyochophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Dinophyceae/Dinophyceae \/ Dinophyceae \/ Dinophyceae/g' tmp #Dinophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Pelagophyceae/Pelagophyceae \/ Pelagophyceae \/ Pelagophyceae/g' tmp #Pelagophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Raphidophyceae/Raphidophyceae \/ Raphidophyceae \/ Raphidophyceae/g' tmp #Raphidophyceae
sed -i 's/unknown kingdom \/ unknown phylum \/ Synurophyceae/Synurophyceae \/ Synurophyceae \/ Synurophyceae/g' tmp #Synurophyceae

#label for those missing kindom + phylum + class labels
sed -i 's/unknown kingdom \/ unknown phylum \/ unknown class \/ Telonemida/Telonemida \/ Telonemida \/ Telonemida \/ Telonemida/g' tmp #Telonemida
sed -i 's/unknown kingdom \/ unknown phylum \/ unknown class \/ Jakobida/Jakobida \/ Jakobida \/ Jakobida \/ Jakobida/g' tmp #Jakobida

#execute final step for the LCA program (forming consensus taxonomies for each ASV)
python2 ~/programs/galaxy-tool-lca/lca.py -i tmp -o blast_96_sim/taxonomy_table.CO1.NCBI_NT+customDB.96sim.txt -b 100 -id 96 -cov 50 -t best_hit -tid 98 -tcov 80 -fh environmental,unidentified,kingdom,phylum,class,order,family,genus -flh unclassified,unknown

#cleanup
mv tmp blast_96_sim/taxonomy_CO1_ASV_sequences.combined.blast.out #replace with taxonomy added blast output
rm taxonomy_tmp blast_96_sim/CO1_ASV_sequences.combined.blast.out blast_96_sim/tmp #remove redundant files


####iterative blast, used with CO1 data from late fall 2021####
#assign taxonomy for dada2-processed CO1 amplicon data with blast using NCBI NT and custom CO1 databases together
mkdir blast_96_sim
# blast against custom blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/CO1_database/blast_DB/CO1.BOLD_genbank_combined.rep_set.blast_DB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.oriented.fasta  -out blast_96_sim/CO1_ASV_sequences.customDB.blast.out
# blast against genbank NT blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2021-11-05/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query CO1_ASV_sequences.oriented.fasta  -out blast_96_sim/CO1_ASV_sequences.blast.out


#filter input fasta to only contain sequences with no hits in the above two blast runs
cat blast_96_sim/CO1_ASV_sequences*blast.out | cut -f1,1 | sort | uniq > blast_96_sim/blast_hit_ASVs.txt
grep -wv -f blast_96_sim/blast_hit_ASVs.txt sequence_ASVname_mapping.txt | cut -f1,1 | sort > blast_96_sim/no_blast_hit_ASVs.txt
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' blast_96_sim/no_blast_hit_ASVs.txt CO1_ASV_sequences.oriented.fasta > blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta

#try blasting this output with lower thresholds for similarity, see what you get
mkdir blast_90_sim
# blast against custom blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 90 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/CO1_database/blast_DB/CO1.BOLD_genbank_combined.rep_set.blast_DB -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta  -out blast_90_sim/CO1_ASV_sequences.customDB.blast.out
# blast against genbank NT blast DB
blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 90 -qcov_hsp_perc 50 -db ~/projects/taxonomyDBs/NCBI_NT/2021-11-05/nt -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta  -out blast_90_sim/CO1_ASV_sequences.blast.out


#IMPORTANT: next steps are done in R, for simplicity's sake. a custom script here would also work. this is simpler.
#we need to add taxonIDs for the customDB (adding them directly to the blast DB has not worked in the past, they don't get returned in the blast output). Using the blast output and a map of accessions to taxonIDs, we add the taxonIDs for each blast result.
library(tidyverse)
library(data.table)
library(ShortRead)
library(Biostrings)
library(seqinr)
taxonidmap <- read.delim("~/projects/taxonomyDBs/CO1_database/taxonID_map/CO1.BOLD_genbank_combined.taxonID_map.w_hakai_barcodes.txt", sep="\t", header=F)
blastfile <- read.delim("blast_96_sim/CO1_ASV_sequences.customDB.blast.out", sep="\t", header=F)
blastfile2 <- read.delim("blast_90_sim/CO1_ASV_sequences.customDB.blast.out", sep="\t", header=F)
blastfile <- rbind(blastfile, blastfile2) #join iterations for customDB blast
colnames(taxonidmap) <- c("accession", "taxonID")
colnames(blastfile) <- c("asv", "col2", "accession", "blasttaxid", "col5", "col6", "col7", "col8")
taxonidmap$accession <- trimws(taxonidmap$accession, which = c("both"))
blastfile_wtaxids <- merge(blastfile,taxonidmap, by="accession", all.x=TRUE)
blastfile_output <- subset(blastfile_wtaxids, select=c("asv", "col2", "accession", "taxonID", "col5", "col6", "col7", "col8")) #it's okay here that "col2" is just all NA values, we need it to conform to the blast output from the NT database, but it doesn't get used by the taxonomy assignment scripts
blastfile_output <- blastfile_output[order(blastfile_output$asv),]
write.table(blastfile_output, "blast_96_sim/CO1_ASV_sequences.customDB_96_90.blast.out", row.names=F, col.names=F, quote=F, sep="\t") #overwriting input

#combine blast results in a way that the add taxonomy and LCA scripts can handle
blastout_customDB <- read.delim("blast_96_sim/CO1_ASV_sequences.customDB_96_90.blast.out", sep="\t", header=F)
blastout_NCBINT <- read.delim("blast_96_sim/CO1_ASV_sequences.blast.out", sep="\t", header=F)
blastout_NCBINT_2 <- read.delim("blast_96_sim/CO1_ASV_sequences.blast.out", sep="\t", header=F)
blastout_NCBINT <- rbind(blastout_NCBINT, blastout_NCBINT_2) #join iterations for NT blast
blastout_combined <- rbind(blastout_customDB, blastout_NCBINT) #combine tables
tmp <- blastout_combined[order(-blastout_combined$V5),] #order descending order for percent identity
tmp <- tmp[order(tmp$V1),] #order by ASV name
blastout_combined <- tmp 
#write to file
write.table(blastout_combined, "blast_96_sim/CO1_ASV_sequences.combined_all.blast.out", row.names=F, col.names=F, quote=F, sep="\t")
#now quit R and continue with the remaining code in your bash shell

#avoid " keyerror: 'NA' " with python script by filtering on "NA" as a wholeword string
#explanation: occasionally, an output line from blast has an NA, which causes an error with the Simple-LCA script below. quick fix is to remove these lines (they're quite rare anyway)
grep -v -w "NA" blast_96_sim/CO1_ASV_sequences.combined_all.blast.out > blast_96_sim/tmp

#execute first step for the LCA program (adding taxonomy strings based on taxonIDs in blast output)
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py -i blast_96_sim/tmp -t ~/projects/taxonomyDBs/NCBI_taxonomy/2021-11-05/rankedlineage.dmp -m ~/projects/taxonomyDBs/NCBI_taxonomy/2021-11-05/merged.dmp -o taxonomy
cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_tmp > tmp

#in this section, a series of taxonomy string modifications are made. This makes the final output more readable/easier to understand, but is optional.
#IMPORTANT: please note that the filtering criteria in the final step depend on some of this filtering (e.g. blast hits with the word "phylum" will be removed. see -fh parameter in final step) 
# because i'm replacing "unknown phylum" in these sed commands, these sequences are retained.
# if you choose not to do the replacement, the blast hits with "unknown plylum" will not be used in LCA determination unless you also change the filtering criteria set in the -fh parameter during the final step.
# also note, "unknown phylum" is present in any taxonomy where the clade's phylum is uncertain in the NCBI taxonomy system, it doesn't indicate any other kind of uncertainty about the provenance of the sequence.

#label fix for clades missing "kingdom" label
sed -i 's/unknown kingdom \/ Bacillariophyta/Eukaryota \/ Bacillariophyta/g' tmp #Bacillariophyta
sed -i 's/unknown kingdom \/ Ciliophora/Eukaryota \/ Ciliophora/g' tmp #Ciliophora
sed -i 's/unknown kingdom \/ Discosea/Eukaryota \/ Discosea/g' tmp #Discosea
sed -i 's/unknown kingdom \/ Evosea/Eukaryota \/ Evosea/g' tmp #Evosea
sed -i 's/unknown kingdom \/ Haptista/Eukaryota \/ Haptista/g' tmp #Haptista
sed -i 's/unknown kingdom \/ Rhodophyta/Eukaryota \/ Rhodophyta/g' tmp #Rhodophyta

#and for those missing kingdom + phylum labels
sed -i 's/Eukaryota \/ unknown phylum \/ Chrysophyceae/Eukaryota \/ Chrysophyceae \/ Chrysophyceae/g' tmp #Chrysophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Cryptophyceae/Eukaryota \/ Cryptophyceae \/ Cryptophyceae/g' tmp #Cryptophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Oomycota/Eukaryota \/ Oomycota \/ Oomycota/g' tmp #Oomycota
sed -i 's/Eukaryota \/ unknown phylum \/ Phaeophyceae/Eukaryota \/ Phaeophyceae \/ Phaeophyceae/g' tmp #Phaeophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Phaeophyceae/Eukaryota \/ Phaeophyceae \/ Phaeophyceae/g' tmp #Phaeophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Bigyra/Eukaryota \/ Bigyra \/ Bigyra/g' tmp #Bigyra
sed -i 's/Eukaryota \/ unknown phylum \/ Dictyochophyceae/Eukaryota \/ Dictyochophyceae \/ Dictyochophyceae/g' tmp #Dictyochophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Dinophyceae/Eukaryota \/ Dinophyceae \/ Dinophyceae/g' tmp #Dinophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Pelagophyceae/Eukaryota \/ Pelagophyceae \/ Pelagophyceae/g' tmp #Pelagophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Raphidophyceae/Eukaryota \/ Raphidophyceae \/ Raphidophyceae/g' tmp #Raphidophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Synurophyceae/Eukaryota \/ Synurophyceae \/ Synurophyceae/g' tmp #Synurophyceae

#label for those missing kindom + phylum + class labels
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Telonemida/Eukaryota \/ Telonemida \/ Telonemida \/ Telonemida/g' tmp #Telonemida
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Jakobida/Eukaryota \/ Jakobida \/ Jakobida \/ Jakobida/g' tmp #Jakobida

#execute final step for the LCA program (forming consensus taxonomies for each ASV)
mv -f tmp blast_96_sim/CO1_ASV_sequences.combined_all.blast.out #put tmp file where it belongs, add label (this is an overwrite, careful!!)
python2 ~/programs/galaxy-tool-lca/lca.py -i blast_96_sim/CO1_ASV_sequences.combined_all.blast.out -o taxonomy_table.CO1.NCBI_NT+customDB.iterative_blast.txt -b 100 -id 80 -cov 50 -t best_hit -tid 98 -tcov 80 -fh environmental,unidentified,kingdom -flh unclassified

#cleanup
rm taxonomy_tmp blast_96_sim/tmp #remove redundant files


