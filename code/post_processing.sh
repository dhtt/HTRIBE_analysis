# TODO: Define path to folder where hypertribe results are stored
wt=( 1 2 3)
mcherry=( 4 5 6)
imp2=( 7 8 9)
declare -a type=("" "_1%" "_5%")


# Compare each replicate of WT group to all samples in mCherry or IMP2 group
# After this step, "wt_comparison" folder contain 18 bedgraph files (3 types * (3 wt replicates  * 2 comparison groups))
OUTDIR=wt_comparison
mkdir -p $OUTDIR/imp2 $OUTDIR/mcherry 
for i in "${wt[@]}"
do
    for j in "${type[@]}"
    do
        bedtools intersect -wa -wb -f 0.9 -r -a IMP2_"$i"_${mcherry[0]}_A2G"$j".bedgraph -b IMP2_"$i"_${mcherry[1]}_A2G"$j".bedgraph IMP2_"$i"_${mcherry[2]}_A2G"$j".bedgraph > $$OUTDIR/wt"$i"_mcherry_A2G"$j".bedgraph
        bedtools intersect -wa -wb -f 0.9 -r -a IMP2_"$i"_${imp2[0]}_A2G"$j".bedgraph -b IMP2_"$i"_${imp2[1]}_A2G"$j".bedgraph IMP2_"$i"_${imp2[2]}_A2G"$j".bedgraph > $$OUTDIR/wt"$i"_imp2_A2G"$j".bedgraph
    done
done


# Subtract WT-mCherry from WT-IMP2
# After this step, "imp2/subtracted" folder contain 9 bedgraph files (3 types * (3 wt replicates)) for IMP2 
cd $OUTDIR
OUTDIR=imp2/subtracted
mkdir -p $OUTDIR
for i in "${wt[@]}"
do
    for j in "${type[@]}"
    do
        bedtools intersect -wa -v -f 0.9 -r -a wt"$i"_imp2_A2G"$j".bedgraph -b wt"$i"_mcherry_A2G"$j".bedgraph > $OUTDIR/HyperTRIBE_wt"$i"_imp2_A2G"$j".bedgraph.temp
        awk 'BEGIN { OFS="\t" } {$25=""; print $0}' $OUTDIR/HyperTRIBE_wt"$i"_imp2_A2G"$j".bedgraph.temp | sed "s/\t\t/\t/g" > $OUTDIR/HyperTRIBE_wt"$i"_imp2_A2G"$j".bedgraph
    done
done



# Combine the results from 3 WT replicates for IMP2 editing sites
# After this step, "imp2/subtracted/combined" folder contain 3 bedgraph files (3 types) for IMP2 
INDIR=$OUTDIR
OUTDIR=imp2/subtracted/combined
mkdir -p $OUTDIR
for i in "${type[@]}"
do
    bedtools intersect -wa -f 0.9 -r -a $INDIR/HyperTRIBE_wt1_imp2_A2G"$i".bedgraph -b $INDIR/HyperTRIBE_wt2_imp2_A2G"$i".bedgraph $INDIR/HyperTRIBE_wt3_imp2_A2G"$i".bedgraph | grep EXON > $OUTDIR/HyperTRIBE_wt_imp2_A2G"$i".bedgraph
done

# Summarize results for IMP2 group 
for HYPERTRIBE_FILE in $OUTDIR/*.bedgraph
do
    OUTFILE="${HYPERTRIBE_FILE%.*}"
    perl $HYPERTRIBE/summarize_results.pl $HYPERTRIBE_FILE > $OUTFILE.xls
done

# Combine the results from 3 WT replicates for mCherry editing sites
# After this step, "imp2/subtracted/combined" folder contain 3 files (3 types) for IMP2 
OUTDIR=mcherry/combined
mkdir -p $OUTDIR
for i in "${type[@]}"
do
    bedtools intersect -wa -f 0.9 -r -a wt1_mcherry_A2G"$i".bedgraph -b wt2_mcherry_A2G"$i".bedgraph wt3_mcherry_A2G"$i".bedgraph > $OUTDIR/HyperTRIBE_wt_mcherry_A2G"$i".bedgraph.temp
    awk 'BEGIN { OFS="\t" } {$25=""; print $0}' $OUTDIR/HyperTRIBE_wt_mcherry_A2G"$i".bedgraph.temp | sed "s/\t\t/\t/g" > $OUTDIR/HyperTRIBE_wt_mcherry_A2G"$i".bedgraph
done

# Summarize results for IMP2 group 
for HYPERTRIBE_FILE in $OUTDIR/*.bedgraph
do
    OUTFILE="${HYPERTRIBE_FILE%.*}"
    perl $HYPERTRIBE/summarize_results.pl $HYPERTRIBE_FILE > $OUTFILE.xls
done

