bedtools subtract -a wt_imp2_A2G_1%.bed -b wt_mcherry_A2G_1%.bed -A > wt_imp2_sub_wt_mcherry_1%.bed
bedtools subtract -a mcherry_imp2_A2G_1%.bed -b wt_mcherry_A2G_1%.bed -A > mcherry_imp2_sub_wt_mcherry_1%.bed
cat mcherry_imp2_sub_wt_mcherry_1%.bed wt_imp2_sub_wt_mcherry_1%.bed | bedtools sort | bedtools merge -c 5,4,5 -o distinct,distinct,count | bedtools sort > subwt_submcherry_1%.bedgraph

cat wt_imp2_A2G_1%.bed mcherry_imp2_A2G_1%.bed | bedtools sort | bedtools merge -c 4,5 -o distinct,distinct > wt.mcherry_imp2_A2G_1%.bed
bedtools subtract -a wt.mcherry_imp2_A2G_1%.bed -b wt_mcherry_A2G_1%.bed -A | bedtools sort | awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$4"\t."}' > wt.mcherry_subimp2_A2G_1%.bedgraph_test

