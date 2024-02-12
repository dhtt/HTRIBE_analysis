from helper_functions import bedtools_sort, bedtools_groupby, bedtools_intersect
from pathlib import Path
import pandas as pd
import os
import plotly.express as px
from plotly.subplots import make_subplots

def gather_meme_motif(STREME_main_path, refgen):
    """
    Gathers meme motifs from multiple comparisons and annotates them with intersecting regions.

    Args:
        STREME_main_path (Path): The main path where the STREME results are stored.
            Default: '~/HTRIBE_analysis/1%_bed_motif/XSTREME/_all/'
        refgen (str): The path to the reference genome file. 
            Default: '~/HTRIBE_analysis/refgen/mm39.ncbiRefSeq.CDS_UTR.gtf'

    Returns:
        None
    """
    comparisons = ['wt_imp2_vs_mcherry_imp2', 'wt_imp2_vs_wt_mcherry', 'mcherry_imp2_vs_wt_mcherry']
    
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


def analyze_meme_motif(STREME_main_path, motives_from_all_comparisons_df='motives_annotated.csv'):
    """
    Analyzes meme motif data and generates a sunburst plot.

    Parameters:
    - STREME_main_path (str): The path to the main directory of the STREME analysis.
    - motives_from_all_comparisons_df (str): The filename of the motives data file. Default is 'motives_annotated.csv'.

    Returns:
    - None
    """
    motives_df = pd.read_csv(STREME_main_path.joinpath(motives_from_all_comparisons_df), header=0, sep='\t')
    
    if os.path.exists(STREME_main_path.joinpath('motives_annotated_freq.csv')):
        motives_df_freq = pd.read_csv(STREME_main_path.joinpath('motives_annotated_freq.csv'), header=0, sep='\t')
    else:
        motives_df_freq = motives_df \
            .groupby(['comparison', 'motif', 'region']).size() \
            .reset_index() 
        motives_df_freq.columns = ['comparison_type', 'motif', 'region', 'count']
        motives_df_freq['motif'] = motives_df_freq['motif'].apply(lambda x: x.split('-')[1])
        motives_df_freq['perc'] = round(100*motives_df_freq['count']/ motives_df_freq.groupby(['comparison_type', 'motif'])['count'].transform('sum'), 2)
        motives_df_freq.to_csv(STREME_main_path.joinpath('motives_annotated_freq.csv'), header=True, sep='\t')
   
    fig = make_subplots(rows=1, cols=3, 
                        specs=[[{'type':'sunburst'}, {'type':'sunburst'}, {'type':'sunburst'}]],
                        subplot_titles=("WT-IMP2<br>vs.<br>WT-mCherry", "mCherry-IMP2<br>vs.<br>WT-mCherry", "WT-IMP2<br>vs.<br>mCherry-IMP2"),
                        horizontal_spacing=0.1)    
    
    plot_colors = {'(?)':'#EDEDED', 'CDS':'#FFD966', '3UTR':'#97D077', '5UTR':'#EA6B66'}
    i = 1
    for part, col in plot_colors.items():
          fig.add_annotation(dict(font=dict(color=col, ),  text='<b>'+part+'</b>', 
                                  y=0, x=1/5*i, showarrow=False, xanchor='center'))
          i+=1
          
          
    fig.add_trace(px.sunburst(motives_df_freq[motives_df_freq['comparison_type'].str.contains('wt_imp2_vs_wt_mcherry')], 
                              values='perc', path=['motif', 'region'], 
                              color='region', color_discrete_map=plot_colors
                        ).data[0], 1, 1)
    
    fig.add_trace(px.sunburst(motives_df_freq[motives_df_freq['comparison_type'].str.contains('mcherry_imp2_vs_wt_mcherry')], 
                              values='perc', path=['motif', 'region'], 
                              color='region', color_discrete_map=plot_colors
                        ).data[0], 1, 2)
    
    fig.add_trace(px.sunburst(motives_df_freq[motives_df_freq['comparison_type'].str.contains('wt_imp2_vs_mcherry_imp2')], 
                              values='perc', path=['motif', 'region'], 
                              color='region', color_discrete_map=plot_colors, 
                        ).data[0], 1, 3)
    
    fig.update_layout_images(xref='paper', yref='paper')
    fig.update_traces(textinfo='percent parent', 
                  #     texttemplate = '%{label}<br>%{percentParent:.1%}',
                      texttemplate = '%{percentParent:.1%}', #_withlabel? True:False
                      insidetextorientation='horizontal',
                      textfont=dict(size=12),
                      marker_line=dict(width=0.5, color='white')
                      )
    fig.update_layout(uniformtext=dict(minsize=9, mode='hide'), #_withlabel? 6:9
                      paper_bgcolor='rgba(0, 0, 0, 0)')
    fig.write_html(STREME_main_path.joinpath('motives_annotated_sunburst_perc.html'))
    fig.write_image(STREME_main_path.joinpath('motives_annotated_sunburst_perc.png'), scale=4)
    
    
if __name__=='__main__':
    import sys
    STREME_main_path = Path(sys.argv[1])
    refgen = sys.argv[2]
    gather_meme_motif(STREME_main_path, refgen)
    analyze_meme_motif(STREME_main_path)
