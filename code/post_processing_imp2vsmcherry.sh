# TODO: Define path to folder where hypertribe results are stored
wt=( 1 2 3)
mcherry=( 4 5 6)
imp2=( 7 8 9)
declare -a type=("" "_1%" "_5%")


# Compare each replicate of WT group to all samples in mCherry or IMP2 group
# After this step, "wt_comparison" folder contain 18 bedgraph files (3 types * (3 wt replicates  * 2 comparison groups))
OUTDIR="imp2_vs_mcherry"
mkdir -p $OUTDIR
for i in "${mcherry[@]}"
do
    for j in "${type[@]}"
    do
        echo "fdsafa"
        bedtools intersect -wa -wb -f 0.9 -r -a IMP2_"$i"_${imp2[0]}_A2G"$j".bedgraph -b IMP2_"$i"_${imp2[1]}_A2G"$j".bedgraph IMP2_"$i"_${imp2[2]}_A2G"$j".bedgraph > $OUTDIR/mcherry"$i"_imp2_A2G"$j".bedgraph
    done
done

# Combine the results from 3 WT replicates for IMP2 editing sites
# After this step, "imp2/subtracted/combined" folder contain 3 bedgraph files (3 types) for IMP2 
# TODO: Define path to folder where hypertribe results are stored
wt=( 1 2 3)
mcherry=( 4 5 6)
imp2=( 7 8 9)
declare -a type=("" "_1%" "_5%")
OUTDIR="imp2_vs_mcherry"
INDIR=$OUTDIR
OUTDIR=imp2_vs_mcherry/combined
mkdir -p $OUTDIR
for i in "${type[@]}"
do
    bedtools intersect -wa -f 0.9 -r -a $INDIR/mcherry4_imp2_A2G"$i".bedgraph -b $INDIR/mcherry5_imp2_A2G"$i".bedgraph $INDIR/mcherry6_imp2_A2G"$i".bedgraph | grep EXON > $OUTDIR/HyperTRIBE_mcherry_imp2_A2G"$i".bedgraph.temp
    awk 'BEGIN { OFS="\t" } {$25=""; print $0}' $OUTDIR/HyperTRIBE_mcherry_imp2_A2G"$i".bedgraph.temp | sed "s/\t\t/\t/g" > $OUTDIR/HyperTRIBE_mcherry_imp2_A2G"$i".bedgraph 
done


declare -a type=("" "_1%" "_5%")
for i in "${type[@]}"
$INDIR='.'
$INDIR='.'
do
    bedtools intersect -u -a ./HyperTRIBE_mcherry_imp2_A2G_5%.bedgraph -b ./HyperTRIBE_wt_imp2_A2G_5%.bedgraph > ./extra/mcherry_imp2_vs_wt_imp2_5%.bedgraph.temp
    bedtools intersect -u -a ./HyperTRIBE_mcherry_imp2_A2G_5%.bedgraph -b ./HyperTRIBE_wt_mcherry_A2G_5%.bedgraph > ./extra/mcherry_imp2_vs_wt_mcherry_5%.bedgraph.temp
    bedtools intersect -u -a ./HyperTRIBE_wt_mcherry_A2G.bedgraph -b ./HyperTRIBE_wt_imp2_A2G_1%.bedgraph > ./extra/wt_mcherry_vs_wt_imp2.bedgraph.temp

    awk 'BEGIN { OFS="\t" } {$25=""; print $0}' ./extra/mcherry_imp2_vs_wt_imp2_5%.bedgraph.temp | sed "s/\t\t/\t/g" > ./extra/mcherry_imp2_vs_wt_imp2_5%.bedgraph
    awk 'BEGIN { OFS="\t" } {$25=""; print $0}' ./extra/mcherry_imp2_vs_wt_mcherry_5%.bedgraph.temp | sed "s/\t\t/\t/g" > ./extra/mcherry_imp2_vs_wt_mcherry_5%.bedgraph
    awk 'BEGIN { OFS="\t" } {$25=""; print $0}' ./extra/wt_mcherry_vs_wt_imp2.bedgraph.temp | sed "s/\t\t/\t/g" > ./extra/wt_mcherry_vs_wt_imp2.bedgraph
done

# Summarize results for IMP2 group 
for HYPERTRIBE_FILE in *.bedgraph
do
    OUTFILE="${HYPERTRIBE_FILE%.*}"
    # bedtools groupby -i $HYPERTRIBE_FILE -g 1-26 -c 1 -o first > $HYPERTRIBE_FILE.grouped
    perl $HYPERTRIBE/summarize_results.pl $HYPERTRIBE_FILE > $OUTFILE.xls
done

rm $OUTDIR/*.temp



wt=( 1 2 3)
mcherry=( 4 5 6)
imp2=( 7 8 9)
declare -a type=("" "_1%" "_5%")
INDIR="."
OUTDIR="./compare_replicates"
ls $INDIR $OUTDIR 
mkdir -p $OUTDIR

for t in "${type[@]}"
do
    for i in "${wt[@]}"  
    # for i in "${mcherry[@]}"
    do
        # for j in "${wt[@]}"
        # for j in "${mcherry[@]}"
        for j in "${imp2[@]}"
        do
            echo $i $j $t
            perl $HYPERTRIBE/summarize_results.pl IMP2_"$i"_"$j"_A2G"$t".bedgraph  > $OUTDIR/test_"$i"_"$j$t".xls
    
        done
    done
done

for ((sample=1;sample<10;sample+=1))
do
    salmon quant -i /home/dhthutrang/TRIBE/refgen/index_salmon_ncbi -l A -1 /home/dhthutrang/TRIBE/mRNA_seq/processed/extract.trim/ctrl_DS"$sample"_R1.fastq -2 /home/dhthutrang/TRIBE/mRNA_seq/processed/extract.trim/ctrl_DS"$sample"_R2.fastq -o alignment_salmon_ncbi/ctrl_DS"$sample" &
done 
