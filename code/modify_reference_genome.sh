
bedtools sort -i /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.transcript.gtf > /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.transcript.gtf.sorted
bedtools merge -i /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.transcript.gtf.sorted -c 3,7,8,9 -o first,first,first,first > /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.transcript.gtf.sorted.merged


bedtools sort -i /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf > /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf.sorted
bedtools merge -i /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf.sorted -c 3,7,8,9 -o distinct,first,first,first > /home/dhthutrang/TRIBE/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf.sorted.merged

