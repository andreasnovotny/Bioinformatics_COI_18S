Rscript Code/HelpScripts/Merge_ASV.r \
    Data/Assign_18S_Taxonomy/ASV \
    Data/18S_QPKbulk_2017

wd='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_18S_Taxonomy'

cd ${wd}

Rscript ../../Code/HelpScripts/AssignTaxonomy_18S.R

cp -r ASV ../../ProcessedData/18S_Andreas
cp tax_tab_18S.RDS ../../ProcessedData/18S_Andreas/taxonomy