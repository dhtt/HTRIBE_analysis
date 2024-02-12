# This script takes two arguments: kmer length and a fasta file.
# It creates a jellyfish count file and generates histogram and count files for each kmer length.

kmer=$1 # Length of kmer from 4 to 8
fa_file=$2 # example: mcherry_imp2.bed.fa
filename=$(basename "$fa_file" .bed.fa) # Remove the suffix .bed.fa from fa_file

jellyfish count -m $kmer -s 100M -t 12 -C $fa_file -o 'jf/'$filename'_'$kmer'.jf'
for file in *$kmer.jf
do 
    jellyfish histo $file -o "${file%.*}"_histo.txt; 
    jellyfish dump $file -o "${file%.*}"_count.txt; 
done;