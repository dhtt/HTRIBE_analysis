background_path=/home/dhthutrang/TRIBE/mRNA_seq/processed/extract.trim.align.dedup/bed

for file in expand_site_nomerged/*/*.bed;
do (
    FILENAME=${file%%.*}
    FILENAME=${FILENAME##*/}
    FILEPATH=${file%/*} 
    echo $FILEPATH/$FILENAME
    mkdir -p $FILEPATH/gene_list
    cut -f5 $file | sort -u > $FILEPATH/gene_list/$FILENAME.txt
    ) 
done    
wait

for file in expand_site_nomerged/*/*.bed;
do (
    FILENAME=${file%%.*}
    FILENAME=${FILENAME##*/}
    FILEPATH=${file%/*} 
    threshold=_${FILENAME##*_}
    if [ $threshold = '_A2G' ]; then
        threshold=''
    fi
    mkdir -p $FILEPATH/findMotifsGenome $FILEPATH/findMotifsGenome_bg $FILEPATH/findMotifsGenome_bg_raw $FILEPATH/findMotifs $FILEPATH/findMotifs_bg $FILEPATH/findMotifs_bg_raw
    
    findMotifsGenome.pl $file mm39 $FILEPATH/findMotifsGenome/$FILENAME/ -bits -len 5 -rna
    findMotifsGenome.pl $file mm39 $FILEPATH/findMotifsGenome_bg/$FILENAME/ -bits -len 5 -bg $FILEPATH/wt_mcherry_A2G$threshold.bed -rna
    findMotifsGenome.pl $file mm39 $FILEPATH/findMotifsGenome_bg_raw/$FILENAME/ -bits -len 5 -bg $background_path/control.bed -rna

    findMotifs.pl $FILEPATH/gene_list/$FILENAME.txt mouse $FILEPATH/findMotifs/$FILENAME -bits -len 5 -rna
    findMotifs.pl $FILEPATH/gene_list/$FILENAME.txt mouse $FILEPATH/findMotifs_bg/$FILENAME -bits -len 5 -bg $FILEPATH/gene_list/wt_mcherry_A2G$threshold.txt -rna
    ) &
done    
wait
echo "FINISH"
