#code for producing taxonomic assignments using blast, for CO1  amplicon data
#author: Evan Morien
#last modified: June 29th, 2021

# Introduction
#	The following code/pipeline was developed over the course of 2020/2021 to use blast to generate taxonomic assignments for both CO1 and 12S amplicon sequencing experiments
#	Code is divided into sections based on the amplicon in question, as well as the reference database

#### CO1 amplicon blast with NT and custom CO1 sequence database ####
#assign taxonomy for dada2-processed CO1 amplicon data with blast using our Hakai CO1 barcode sequences



# Start by merging ASV tables from DADA2
# CHANGE ME: to output directory followed by several DADA analysis directories:
Rscript Code/HelpScripts/Merge_ASV.r \
    Data/Assign_COI_Taxonomy/ASV \
    Data/COI_QU39-2017 \
    Data/COI_QU39-2018 \
    Data/COI_QU39-2019 \
    Data/COI_Zoopsprint2022 \
    Data/COI_QPKbulk_2017
    

# CHANGE ME: Define paths to databases, querry sequences (produced by script above) and 
BOLD_genbank_combo='/data/taxonomyDBs/CO1_database/blast_DB/CO1.BOLD_genbank_combined.rep_set.blast_DB'
genbank_NT_blast='/data/taxonomyDBs/NCBI_NT/2023-02-07/nt'
querry='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_COI_Taxonomy/ASV'
wd='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_COI_Taxonomy'


#############################################################

cd ${wd}

####iterative blast, used with CO1 data from late fall 2021####
#assign taxonomy for dada2-processed CO1 amplicon data with blast using NCBI NT and custom CO1 databases together
mkdir blast_96_sim

# blast against custom blast DB
blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 96 \
    -qcov_hsp_perc 50 \
    -db ${BOLD_genbank_combo} \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query ${querry}/CO1_ASV_sequences.fasta  \
    -out blast_96_sim/CO1_ASV_sequences.customDB.blast.out


# blast against genbank NT blast DB
blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 96 \
    -qcov_hsp_perc 50 \
    -db ${genbank_NT_blast} \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query ${querry}/CO1_ASV_sequences.fasta  \
    -out blast_96_sim/CO1_ASV_sequences.blast.out



#filter input fasta to only contain sequences with no hits in the above two blast runs
cat blast_96_sim/CO1_ASV_sequences*blast.out | cut -f1,1 | sort | uniq > blast_96_sim/blast_hit_ASVs.txt
grep -wv -f blast_96_sim/blast_hit_ASVs.txt ${querry}/sequence_ASVname_mapping.CO1.txt | cut -f1,1 | sort > blast_96_sim/no_blast_hit_ASVs.txt
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' blast_96_sim/no_blast_hit_ASVs.txt ${querry}/CO1_ASV_sequences.fasta > blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta

#blast this output with lower thresholds for similarity
mkdir blast_90_sim

# blast against custom blast DB
blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 90 \
    -qcov_hsp_perc 50 \
    -db ${BOLD_genbank_combo} \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta  \
    -out blast_90_sim/CO1_ASV_sequences.customDB.blast.out

# blast against genbank NT blast DB
blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 90 \
    -qcov_hsp_perc 50 \
    -db ${genbank_NT_blast} \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query blast_96_sim/CO1_ASV_sequences.no_blast_hit.fasta  \
    -out blast_90_sim/CO1_ASV_sequences.blast.out


#filter input fasta to only contain sequences with no hits in the above two blast runs
cat blast_90_sim/CO1_ASV_sequences*blast.out | cut -f1,1 | sort | uniq > blast_90_sim/blast_hit_ASVs.txt
cat blast_90_sim/blast_hit_ASVs.txt blast_96_sim/blast_hit_ASVs.txt > blast_90_sim/blast_90_96_hit_ASVs.txt
grep -wv -f blast_90_sim/blast_90_96_hit_ASVs.txt ${querry}/sequence_ASVname_mapping.CO1.txt | cut -f1,1 | sort > blast_90_sim/no_blast_hit_90_96_ASVs.txt
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' blast_90_sim/no_blast_hit_90_96_ASVs.txt ${querry}/CO1_ASV_sequences.fasta > blast_90_sim/CO1_ASV_sequences.no_blast_hit.fasta


#blast this output with lower thresholds for similarity
mkdir blast_80_sim
# blast against custom blast DB
blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 80 \
    -qcov_hsp_perc 50 \
    -db ${BOLD_genbank_combo} \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query blast_90_sim/CO1_ASV_sequences.no_blast_hit.fasta \
    -out blast_80_sim/CO1_ASV_sequences.customDB.blast.out

# blast against genbank NT blast DB
blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 80 \
    -qcov_hsp_perc 50 \
    -db ${genbank_NT_blast} \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query blast_90_sim/CO1_ASV_sequences.no_blast_hit.fasta  \
    -out blast_80_sim/CO1_ASV_sequences.blast.out

cat blast_80_sim/CO1_ASV_sequences*blast.out | cut -f1,1 | sort | uniq > blast_80_sim/blast_hit_ASVs.txt
cat blast_80_sim/blast_hit_ASVs.txt blast_90_sim/blast_hit_ASVs.txt blast_96_sim/blast_hit_ASVs.txt > blast_80_sim/blast_80_90_96_hit_ASVs.txt
grep -wv -f blast_80_sim/blast_80_90_96_hit_ASVs.txt ${querry}/sequence_ASVname_mapping.CO1.txt | cut -f1,1 | sort > blast_80_sim/no_blast_hit_80_90_96_ASVs.txt
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' blast_80_sim/no_blast_hit_80_90_96_ASVs.txt ${querry}/CO1_ASV_sequences.fasta > blast_90_sim/CO1_ASV_sequences.no_blast_hit.fasta

#IMPORTANT: next steps are done in R,
#we need to add taxonIDs for the customDB (adding them directly to the blast DB has not worked in the past, they don't get returned in the blast output). Using the blast output and a map of accessions to taxonIDs, we add the taxonIDs for each blast result.
Rscript ../../Code/HelpScripts/modify.blast.COI.R

#avoid " keyerror: 'NA' " with python script by filtering on "NA" as a wholeword string
#explanation: occasionally, an output line from blast has an NA, which causes an error with the Simple-LCA script below. quick fix is to remove these lines (they're quite rare anyway)
grep -v -w "NA" blast_96_sim/CO1_ASV_sequences.combined_all.blast.out > blast_96_sim/tmp

#execute first step for the LCA program (adding taxonomy strings based on taxonIDs in blast output)
python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py \
    -i blast_96_sim/tmp \
    -t /data/taxonomyDBs/NCBI_taxonomy/2023-02-03/rankedlineage.dmp \
    -m /data/taxonomyDBs/NCBI_taxonomy/2023-02-03/merged.dmp \
    -o taxonomy #HERE

cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) taxonomy_tmp > tmp

#in this section, a series of taxonomy string modifications are made. This makes the final output more readable/easier to understand, but is optional.
#IMPORTANT: please note that the filtering criteria in the final step depend on some of this filtering (e.g. blast hits with the word "phylum" will be removed. see -fh parameter in final step) 
# because i'm replacing "unknown phylum" in these sed commands, these sequences are retained.
# if you choose not to do the replacement, the blast hits with "unknown plylum" will not be used in LCA determination unless you also change the filtering criteria set in the -fh parameter during the final step.
# also note, "unknown phylum" is present in any taxonomy where the clade's phylum is uncertain in the NCBI taxonomy system, it doesn't indicate any other kind of uncertainty about the provenance of the sequence

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
sed -i 's/Eukaryota \/ unknown phylum \/ Bolidophyceae/Eukaryota \/ Bolidophyceae \/ Bolidophyceae/g' tmp #Bolidophyceae
sed -i 's/Eukaryota \/ unknown phylum \/ Polycystinea/Eukaryota \/ Polycystinea \/ Polycystinea/g' tmp #Polycystinea
sed -i 's/Eukaryota \/ unknown phylum \/ Choanoflagellata/Eukaryota \/ Choanoflagellata \/ Choanoflagellata/g' tmp #Choanoflagellata
sed -i 's/Eukaryota \/ unknown phylum \/ Filasterea/Eukaryota \/ Filasterea \/ Filasterea/g' tmp #Filasterea

#label for those missing kindom + phylum + class labels
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Telonemida/Eukaryota \/ Telonemida \/ Telonemida \/ Telonemida/g' tmp #Telonemida
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Jakobida/Eukaryota \/ Jakobida \/ Jakobida \/ Jakobida/g' tmp #Jakobida
sed -i 's/Eukaryota \/ unknown phylum \/ unknown class \/ Pirsoniales/Eukaryota \/ Pirsoniales \/ Pirsoniales \/ Pirsoniales/g' tmp #Pirsoniales

#execute final step for the LCA program (forming consensus taxonomies for each ASV)
mv -f tmp blast_96_sim/CO1_ASV_sequences.combined_all.blast.out #put tmp file where it belongs, add label (this is an overwrite, careful!!)


# python2 ~/programs/galaxy-tool-lca/lca.py \
#     -i blast_96_sim/CO1_ASV_sequences.combined_all.blast.out \
#     -o taxonomy_table.iterative_blast.txt \
#     -b 100 \
#     -id 80 \
#     -cov 50 \
#     -t best_hit \
#     -tid 98 \
#     -tcov 50 \
#     -fh environmental,unidentified,kingdom \
#     -flh unclassified

# #cleanup
# rm taxonomy_tmp blast_96_sim/tmp #remove redundant files

# lca.species.py
# -i	Input file
# -o	Output file
# -b	Bitscore top percentage threshold
# -id	Minimum identity threshold
# -cov	Minimum coverage threshold
# -t	Check the top hit first or perform an lca analysis on all hits. Options:['only_lca', 'best_hit', "best_hits_range"]
# -tid	Identity threshold for the tophit, only used when -t best_hit or best_hits_range
# -tcov	Coverage threshold for the tophit, only used when -t best_hit or best_hits_range
# -fh	Filter hits, filter out lines that contain a certain string
# -flh	Filter lca hits, during the determination of the lca ignore this string
# -minbit	Minimum bitscore threshold


 python2 ~/programs/galaxy-tool-lca/lca.species.py \
     -i blast_96_sim/CO1_ASV_sequences.combined_all.blast.out \
     -o taxonomy_table.CO1.iterative_blast.LCA+best_hit.txt \
     -b 100 \
     -id 90 \
     -cov 50 \
     -t best_hit \
     -tid 98 \
     -tcov 50 \
     -fh environmental,unidentified,kingdom \
     -flh unclassified

 python2 ~/programs/galaxy-tool-lca/lca.species.py \
     -i blast_96_sim/CO1_ASV_sequences.combined_all.blast.out \
     -o taxonomy_table.CO1.iterative_blast.LCA_only.txt \
     -b 100 \
     -id 90 \
     -cov 50 \
     -t only_lca \
     -fh environmental,unidentified,kingdom \
     -flh unclassified




#Move final outputs away from gitignore space
cp taxonomy_table.CO1.iterative_blast.LCA+best_hit.txt ../../ProcessedData/COI_Andreas/COI_taxonomy
cp taxonomy_table.CO1.iterative_blast.LCA_only.txt ../../ProcessedData/COI_Andreas/COI_taxonomy
cp blast_96_sim/CO1_ASV_sequences.combined_all.blast.out ../../ProcessedData/COI_Andreas/COI_taxonomy

cp -r ASV ../../ProcessedData/COI_Andreas