Rscript Code/HelpScripts/Merge_ASV.r \
    Data/Assign_18S_Taxonomy/ASV \
    Data/18S_QPKbulk_2017 \
    Data/18S_QU39_1 \
    Data/18S_QU39_2 \
    Data/18S_QU39_3 \
    Data/18S_QU39_4 \
    Data/18S_QU39_5

wd='/home/andreas.novotny/AmpliconSeqAnalysis/Data/Assign_18S_Taxonomy'

cd ${wd}

Rscript ../../Code/HelpScripts/AssignTaxonomy_18S.R

cp -r ASV ../../ProcessedData/18S_Andreas
cp tax_tab_18S_pr2.RDS ../../ProcessedData/18S_Andreas/tax_tab_18S_pr2.RDS
cp tax_tab_18S_SILVA.RDS ../../ProcessedData/18S_Andreas/tax_tab_18S_SILVA.RDS