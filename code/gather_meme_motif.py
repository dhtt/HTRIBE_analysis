#conda env python3.11
'''
This script gather all MEME-suite/XSTREME motives and annotate them to regions from a reference genome gtf file.

Usage: gather_meme_motif.py 
'''

from helper_functions import bedtools_sort, bedtools_groupby, bedtools_intersect
from pathlib import Path
import pandas as pd

if __name__=='__main__':
    comparisons = ['wt_imp2_vs_mcherry_imp2', 'wt_imp2_vs_wt_mcherry', 'mcherry_imp2_vs_wt_mcherry']
    STREME_main_path = Path('/home/dhthutrang/HTRIBE_analysis/1%_bed_motif/XSTREME/_all/')
    refgen = '/home/dhthutrang/HTRIBE_analysis/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf'
    
    motives_from_all_comparisons = []
    for comparison in comparisons:
        STREME_result_path = STREME_main_path.joinpath(comparison, 'streme_out')
        
        sites_filename = 'sites.tsv'
        sites_df = pd.read_csv(STREME_result_path.joinpath(sites_filename), header=0, sep='\t') \
            .dropna(axis=0, how='any') \
            .astype({'site_Start': int, 'site_End': int}) \
            .groupby(by='motif_ID', axis=0)
        motives = sites_df.groups.keys()
        
        for i, motif in enumerate(motives):
            motif_bed_file = str(STREME_result_path.joinpath(motif + '.bed'))
            motif_bed_file_annotated = motif_bed_file + '.annot'
            motif_df = sites_df \
                .get_group(motif) \
                .drop(['motif_ID', 'motif_ALT_ID'], axis=1) \
                .to_csv(motif_bed_file, index=False, header=False, sep='\t')
        
            bedtools_intersect(bed_file_1=motif_bed_file,
                               bed_file_2=refgen,
                               output_path=motif_bed_file_annotated,
                               intersect_extra_args='-wa -wb')
            
            bedtools_sort(bed_file=motif_bed_file_annotated, 
                          extra_sort_args=' -k5,5n ')
            
            bedtools_groupby(bed_file=motif_bed_file_annotated, 
                             group_by='1,2,3,4,9', 
                             retained_columns='9', 
                             group_mode='distinct')
            
            annotated_regions = pd.read_csv(motif_bed_file_annotated, header=None, sep='\t',
                                            names=['chr', 'start', 'end', 'strand', 'region', 'region_'])
            annotated_regions['motif'] = motif
            annotated_regions['comparison_type'] = comparison
            annotated_regions = annotated_regions.drop(['region_'], axis=1)
            
            motives_from_all_comparisons.append(annotated_regions)
    
    motives_from_all_comparisons_df = pd.concat(motives_from_all_comparisons)
    motives_from_all_comparisons_df.to_csv(STREME_main_path.joinpath('motives_annotated.csv'), header=True, sep='\t')