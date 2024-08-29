SPAN=(100 250 500)
for span in ${SPAN[@]}; do
    echo "Span: "$span
    mkdir -p expand_site_merged/$span
    mkdir -p expand_site_nomerged/$span
    for file in *.bedgraph.site; do
        bedtools merge -c 4,8 -o max,distinct -i $file | bedtools slop -b $span -g /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.dexseq.gtf | bedtools merge -c 4,5 -o max,distinct > expand_site_merged/$span/${file%%.*}.bed
        bedtools merge -c 4,8 -o max,distinct -i $file | bedtools slop -b $span -g /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.dexseq.gtf > expand_site_nomerged/$span/${file%%.*}.bed
    done
done
