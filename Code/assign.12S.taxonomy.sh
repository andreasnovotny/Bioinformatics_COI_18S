
# Start by merging ASV tables from DADA2
# CHANGE ME: to output directory followed by several DADA analysis directories:
Rscript Code/HelpScripts/Merge_ASV.r \
    Data/Assign_12S_Taxonomy/ASV \
    Data/12S_QU39_2018

# CHANGE ME: Define paths to databases, querry sequences (produced by script above) and 

genbank_NT_blast='/data/taxonomyDBs/NCBI_NT/2023-02-07/nt'
querry='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_12S_Taxonomy/ASV'
wd='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_12S_Taxonomy'

#########################

cd ${wd}

#### 12S or 16S blasting against NCBI NT ####
#remember to change the input and output file names so they match the amplicon you're working with
#assign taxonomy with blast NT database at 96% similarity threshold using both 'LCA + besthit' and 'LCA only' parameters
mkdir blast_96_sim_LCA_besthit

blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 96 \
    -qcov_hsp_perc 50 \
    -db ${genbank_NT_blast} \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query ${querry}/CO1_ASV_sequences.fasta  \
    -out blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out


python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py \
    -i blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out \
    -t /data/taxonomyDBs/NCBI_taxonomy/2023-02-03/rankedlineage.dmp \
    -m /data/taxonomyDBs/NCBI_taxonomy/2023-02-03/merged.dmp \
    -o blast_96_sim_LCA_besthit/taxonomy

cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) \
    taxonomy_12S_ASV_sequences.length_var.blast.out > tmp

python2 ~/programs/galaxy-tool-lca/lca.species.py \
    -i tmp \
    -o blast_96_sim_LCA_besthit/taxonomy_table.12S.NCBI_NT.96sim.LCA+besthit.txt \
    -b 100 \
    -id 96 \
    -cov 50 \
    -t best_hit \
    -tid 98 \
    -tcov 80 \
    -fh environmental,unidentified,kingdom \
    -flh unclassified

mkdir blast_96_sim_LCA_only

python2 ~/programs/galaxy-tool-lca/lca.species.py \
    -i tmp \
    -o blast_96_sim_LCA_only/taxonomy_table.12S.NCBI_NT.96sim.LCA_only.txt \
    -b 100 \
    -id 96 \
    -cov 50 \
    -t only_lca \
    -fh environmental,unidentified \
    -flh unclassified

#cleanup
rm blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out #remove blast output without taxonomy
rm taxonomy_12S_ASV_sequences.length_var.blast.out #remove redundant file
mv tmp blast_96_sim_LCA_besthit/12S_ASV_sequences.length_var.blast.out #replace with taxonomy added blast output











# CHANGE ME: Define paths to databases, querry sequences (produced by script above) and 
NCBI_SSU_euks='/data/taxonomyDBs/NCBI_NT/2023-02-07/nt'
querry='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_18S_Taxonomy/ASV'
wd='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_18S_Taxonomy'

cd ${wd}


########
NCBI_SSU_euks='/data/taxonomyDBs/NCBI_SSU_euks/SSU_eukaryote_rRNA '
genbank_NT_blast='/data/taxonomyDBs/NCBI_NT/2023-02-07/nt'
BOLD_genbank_combo='/data/taxonomyDBs/CO1_database/blast_DB/CO1.BOLD_genbank_combined.rep_set.blast_DB'
querry='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_18S_Taxonomy/ASV'
wd='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_18S_Taxonomy'
cd ${wd}
mkdir blast_96_sim

# blast against custom blast DB
blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 96 \
    -qcov_hsp_perc 50 \
    -db ${genbank_NT_blast} \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query ${querry}/CO1_ASV_sequences.fasta  \
    -out blast_96_sim/CO1_ASV_sequences.customDB.blast.out
#########


#### 18S amplicon blast with NCBI SSU database ####
mkdir blast_96_sim

blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 96 \
    -qcov_hsp_perc 50 \
    -db /data/taxonomyDBs/NCBI_SSU_euks/SSU_eukaryote_rRNA \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query ASV/CO1_ASV_sequences.fasta  \
    -out blast_96_sim/18S_ASV_sequences.NCBI_SSU.blast.out


blastn -task megablast -num_threads 38 -evalue 1e-5 -max_target_seqs 10 -perc_identity 96 -qcov_hsp_perc 50 -db /data/taxonomyDBs/NCBI_SSU_euks/SSU_eukaryote_rRNA -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' -query ASV/CO1_ASV_sequences.fasta  -out blast_96_sim/18S_ASV_sequences.NCBI_SSU.blast.out

blastn -task megablast \
    -num_threads 38 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -perc_identity 96 \
    -qcov_hsp_perc 50 \
    -db /data/taxonomyDBs/NCBI_SSU_euks/SSU_eukaryote_rRNA \
    -outfmt '6 qseqid stitle sacc staxid pident qcovs evalue bitscore' \
    -query ASV/CO1_ASV_sequences.fasta  \
    -out blast_96_sim/CO1_ASV_sequences.customDB.blast.out

python2 ~/programs/galaxy-tool-BLAST/blastn_add_taxonomy_lite.py \ #ok
    -i blast_96_sim/18S_ASV_sequences.NCBI_SSU.blast.out \ #ok
    -t /data/taxonomyDBs/NCBI_taxonomy/2023-02-03/rankedlineage.dmp \
    -m /data/taxonomyDBs/NCBI_taxonomy/2023-02-03/merged.dmp -o taxonomy #HERE

cat <(head -n 1 ~/programs/galaxy-tool-lca/example/example.tabular) \ #ok
    taxonomy_18S_ASV_sequences.NCBI_SSU.blast.out > tmp #FIX

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

#for 18S, i have found that the filtering parameters below are best for identifying bacterial (non 18S sequence)
python2 ~/programs/galaxy-tool-lca/lca.species.py \
    -i tmp -o taxonomy_table.18S.NCBI_SSU.96sim.txt \
    -b 100 -id 96 \
    -cov 50 \
    -t best_hit \
    -tid 98 \
    -tcov 80 \
    -fh environmental,unidentified,kingdom \
    -flh unclassified

#cleanup
rm taxonomy_18S_ASV_sequences.NCBI_SSU.blast.out tmp