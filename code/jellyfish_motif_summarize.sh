jellyfish count -m 5 -s 100M -t 12 -C mcherry_imp2.bed.fa -o jf/mcherry_imp2_5.jf
for file in *4.jf
do 
    filename="${file%.*}"; jellyfish histo $file -o ${filename}_histo.txt; 
    filename="${file%.*}"; jellyfish dump $file -o ${filename}_count.txt; 
done;