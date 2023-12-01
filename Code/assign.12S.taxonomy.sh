
# Start by merging ASV tables from DADA2
# CHANGE ME: to output directory followed by several DADA analysis directories:
Rscript Code/HelpScripts/Merge_ASV.r \
    Data/Assign_12S_Taxonomy/ASV \
    Data/12S_QU39_2018 \
    Data/12S_QU39_2017 \
    Data/12S_QU39_2016

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

#Move final outputs away from gitignore space
cp blast_96_sim_LCA_besthit/taxonomy_table.12S.NCBI_NT.96sim.LCA+besthit.txt ../../ProcessedData/12S_Andreas/Taxonomy
cp blast_96_sim_LCA_only/taxonomy_table.12S.NCBI_NT.96sim.LCA_only.txt ../../ProcessedData/12S_Andreas/Taxonomy
cp -r ASV ../../ProcessedData/12S_Andreas



