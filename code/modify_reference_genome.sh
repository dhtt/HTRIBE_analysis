# This script performs sorting and merging operations on two GTF files using bedtools commands subsequently.

bedtools sort -i ~/TRIBE/refgen/mm39.ncbiRefSeq.transcript.gtf > ~/TRIBE/refgen/mm39.ncbiRefSeq.transcript.gtf.sorted
bedtools merge -i ~/TRIBE/refgen/mm39.ncbiRefSeq.transcript.gtf.sorted -c 3,7,8,9 -o first,first,first,first > ~/TRIBE/refgen/mm39.ncbiRefSeq.transcript.gtf.sorted.merged

bedtools sort -i ~/TRIBE/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf > ~/TRIBE/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf.sorted
bedtools merge -i ~/TRIBE/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf.sorted -c 3,7,8,9 -o distinct,first,first,first > ~/TRIBE/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf.sorted.merged
