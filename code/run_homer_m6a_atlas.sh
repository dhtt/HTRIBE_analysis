# m6A-Atlas2_Mouse_exomePeak2.txt and m6A-Atlas2_Mouse_Site.txt are from m6A-Atlas2 database

workdir=~/TRIBE/mRNA_seq/processed/extract.trim.align.dedup/test1/result/diff_span/all_span/CDS_UTR/OR/OR/expand_site_nomerged/m6a_atlas
cd $workdir

awk -v OFS='\t' '{print $2, $3, $4, $7, $11}' m6A-Atlas2_Mouse_exomePeak2.txt > m6aexome.bed
awk -v FS='\t' -v OFS='\t' '{print $2,$3, $3, $6}' m6A-Atlas2_Mouse_Site.txt > m6asite.bed

cat *_A2G_1%.bed > all_A2G_1%.bed

# Reversely, look for motifs in the HTRIBE sites
bedtools intersect -a all_A2G_1%.bed -b m6aexome.bed -wa -wb > common_m6a.bed && sortBed -i common_m6a.bed > common_m6a.sorted.bed && bedtools merge -i common_m6a.sorted.bed -c 1,2,3,4,5,7,8,9,10 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct > common_m6a.sorted.merged.bed
findMotifsGenome.pl common_m6a.bed mm39 m6a_atlas_exome_common/ -bits -len 5 -rna

bedtools intersect -a all_A2G_1%.bed -b m6asite.bed -wa -wb > common_m6a.bed && sortBed -i common_m6a.bed > common_m6a.sorted.bed && bedtools merge -i common_m6a.sorted.bed -c 1,2,3,4,5,7,8,9,10 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct,distinct > common_m6a.sorted.merged.bed
findMotifsGenome.pl common_m6a.bed mm39 m6a_atlas_common/ -bits -len 5 -rna

# Reversely, look for motifs in the exome peaks
bedtools intersect -a all_A2G_1%.bed -b m6aexome.bed -wb > common_m6a.bed && sortBed -i common_m6a.bed > common_m6a.sorted.bed && bedtools merge -i common_m6a.sorted.bed -c 1,2,3,4,5,9,10 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct > common_m6a.sorted.merged.bed
findMotifsGenome.pl common_m6a.bed mm10 m6a_atlas_exome_common/ -bits -len 5 -rna

bedtools intersect -a all_A2G_1%.bed -b m6asite.bed -wb > common_m6a.bed && sortBed -i common_m6a.bed > common_m6a.sorted.bed && bedtools merge -i common_m6a.sorted.bed -c 1,2,3,4,5,9,10 -o distinct,distinct,distinct,distinct,distinct,distinct,distinct > common_m6a.sorted.merged.bed
findMotifsGenome.pl common_m6a.bed mm10 m6a_atlas_common/ -bits -len 5 -rna
