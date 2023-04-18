#!/bin/sh

HyperTRIBE_DIR="~/TRIBE/HTRIBE_analysis/code/HyperTRIBE"

# The combination of tablename, expt name and tp value is used to extract the base composition at a given location of the genome for each experimental condition. the find_rnaeditsites.pl script uses these variables to extract and compare the base composition between the rna library and gDNA library to call the edit sites
#---------------------------
# edit the following varibales as need
annotationfile=$HyperTRIBE_DIR"/annotations/refFlat.txt"

wtRNAtablename=" RNA2 "
wtRNAexp="IMP2"
wtRNAtp="1"
RNAtablename="RNA2"
RNAexp="IMP2"

#edit the timepoint array as needed
# timepoint=(3 4 5)
timepoint=(2 3)
# one or more samples can be run by altering this array. 
#---------------

for tp in ${timepoint[@]}
do
  outfile1=$RNAexp"_"$wtRNAtp"_"$tp"_A2G.txt"
  perl $HyperTRIBE_DIR/find_rnaeditsites.pl -a $annotationfile -t $RNAtablename -e $RNAexp -c $tp -o $outfile1 -g $wtRNAtablename -j $wtRNAexp -k $wtRNAtp  
  
  # Previously used Threshold_editsites_20reads.py script is replaced with filter_by_threshold_without_header.pl for usability reasons and we are using a lower threshold of 5% for editing because empirical data from our lab has shown that the  target list does not change much the threshold is lower. 

  # convert to bedgraph format
  prefix=${outfile1%.txt*}
  python $HyperTRIBE_DIR/convert_editsites_to_bedgraph.py $outfile1
  mv $outfile1".bedgraph" $prefix".bedgraph"

  # apply edit % and read threshold, the index is zero based in the perl script
  # 4th col is edit thresold: $threshold
  # 21st column has a read threshold of 10
  # create a 5% threshold edit file
  edit_threshold=5
  read_threshold=20
  prefix=${outfile1%.txt*}
  out_bedgraph=$prefix"_"$edit_threshold"%.bedgraph"
  perl $HyperTRIBE_DIR/filter_by_threshold_without_header.pl 3 $edit_threshold 20 $read_threshold $prefix".bedgraph" > $out_bedgraph

  # create a 1% threshold file
  edit_threshold=1
  read_threshold=20
  prefix=${outfile1%.txt*}
  out_bedgraph=$prefix"_"$edit_threshold"%.bedgraph"
  perl $HyperTRIBE_DIR/filter_by_threshold_without_header.pl 3 $edit_threshold 20 $read_threshold $prefix".bedgraph" > $out_bedgraph

  #Optional create the target list without removing the background 
  #perl $HyperTRIBE_DIR/summarize_results.pl $out_bedgraph > $prefix"_"$edit_threshold"%_results.xls"

# ----------------------------


done
