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
    mkdir -p $FILEPATH/findMotifsGenome $FILEPATH/findMotifsGenome_bg $FILEPATH/findMotifs $FILEPATH/findMotifs_bg
    
    findMotifsGenome.pl $file mm39 $FILEPATH/findMotifsGenome/$FILENAME/ -bits -len 5
    findMotifsGenome.pl $file mm39 $FILEPATH/findMotifsGenome_bg/$FILENAME/ -bits -len 5 -bg $FILEPATH/wt_mcherry_A2G_$threshold%.bed
    findMotifs.pl gene_list/$FILENAME.txt mouse $FILEPATH/findMotifs/$FILENAME -bits -len 5
    echo "findMotifs.pl $FILEPATH/gene_list/$FILENAME.txt mouse $FILEPATH/findMotifs_bg/$FILENAME -bits -len 5 -bg $FILEPATH/gene_list/wt_mcherry_A2G$threshold.txt"
    ) &
done    

