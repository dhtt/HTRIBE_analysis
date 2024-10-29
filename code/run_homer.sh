background_path=~/TRIBE/mRNA_seq/processed/extract.trim.align.dedup/bed
bed_path=~/TRIBE/mRNA_seq/processed/extract.trim.align.dedup/test1/result/diff_span/all_span/CDS_UTR/OR/OR/expand_site_merged
cd $bed_path
for file in 500/*.bed;
do (
    FILENAME=${file%%.*}
    FILENAME=${FILENAME##*/}
    FILEPATH=${file%/*} 
    mkdir -p $FILEPATH/gene_list
    cut -f5 $file | sort -u > $FILEPATH/gene_list/$FILENAME.txt
    ) 
done    
wait

for file in 500/*.bed;
do (
    FILENAME=${file%%.*}
    FILENAME=${FILENAME##*/}
    FILEPATH=${file%/*} 
    threshold=_${FILENAME##*_}
    if [ $threshold = '_A2G' ]; then
        threshold=''
    fi
    echo $FILEPATH/$FILENAME
    mkdir -p $FILEPATH/findMotifsGenome $FILEPATH/findMotifsGenome_bg $FILEPATH/findMotifsGenome_bg_raw $FILEPATH/findMotifs $FILEPATH/findMotifs_bg $FILEPATH/findMotifs_bg_raw
    
    grep 3UTR $file > $file.3UTR
    grep 5UTR $file > $file.5UTR
    grep CDS $file > $file.CDS

    REG=(3UTR 5UTR CDS)
    for region in ${REG[@]}; 
    do 
        findMotifsGenome.pl $file.$region mm39 $FILEPATH/findMotifsGenome/$FILENAME.$region/ -bits -len 5 -rna
        findMotifsGenome.pl $file.$region mm39 $FILEPATH/findMotifsGenome_bg/$FILENAME.$region/ -bits -len 5 -bg $FILEPATH/wt_mcherry_A2G$threshold.bed -rna
    done

    findMotifsGenome.pl $file mm39 $FILEPATH/findMotifsGenome/$FILENAME/ -bits -len 5 -rna
    findMotifsGenome.pl $file mm39 $FILEPATH/findMotifsGenome_bg/$FILENAME/ -bits -len 5 -bg $FILEPATH/wt_mcherry_A2G$threshold.bed -rna

    findMotifsGenome.pl $file mm39 $FILEPATH/findMotifsGenome_bg_raw/$FILENAME/ -bits -len 5 -bg $background_path/control.bed -rna

    findMotifs.pl $FILEPATH/gene_list/$FILENAME.txt mouse $FILEPATH/findMotifs/$FILENAME -bits -len 5 -rna
    findMotifs.pl $FILEPATH/gene_list/$FILENAME.txt mouse $FILEPATH/findMotifs_bg/$FILENAME -bits -len 5 -bg $FILEPATH/gene_list/wt_mcherry_A2G$threshold.txt -rna
    ) &
done    
wait
echo "FINISH"
