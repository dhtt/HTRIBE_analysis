#!/usr/bin/env python # [1]
'''
This script analyzes the HYPERTRIBE summarized results for each threshold generated (0%, 1%, 5%)

Usage: analyze_by_threshold.py HTRIBE_result_path threshold
'''

from helper_functions import read_bedgraph, create_new_folder
from pathlib import Path
import pandas as pd
from sklearn.decomposition import PCA
import plotly.express as px
import plotly.graph_objects as go 
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
from itertools import combinations

if __name__=='__main__':
     # HTRIBE_result_path = sys.argv[1]
    # threshold = sys.argv[2]

    HTRIBE_result_path = Path('/home/dhthutrang/TRIBE/mRNA_seq/processed/extract.trim.align.dedup/test1/result/collapsed_replicates/comparisons/combined_result/')
    HTRIBE_annot_result_path = HTRIBE_result_path.joinpath('annotated')
    threshold = '5%'
    comparison_type = ['wt_mcherry', 'wt_imp2', 'mcherry_imp2']
    comparison_pair = list(combinations(comparison_type, 2))
    
    # Create a folder storing results
    analyzed_result_path = Path(HTRIBE_result_path).joinpath('analyzed_result')
    create_new_folder(analyzed_result_path)

    result_df = read_bedgraph(result_dir=HTRIBE_annot_result_path, threshold=threshold, file_extension='.bed', header=None)
    result_df.columns = ['chr', 'start', 'end', 'perc', 'id', 'region', 'comparison_type']
    result_df['region'] = result_df['region'].str.split(pat=',')
    result_df = result_df.explode('region')
    result_df_count = result_df.groupby(['comparison_type', 'region'], as_index=False).count().iloc[:, :3]
    result_df_count.columns = ['comparison_type', 'region', 'count']
    result_df_count['perc'] = round(result_df_count['count']/sum(result_df_count['count'])*100, 2)
    
    
    fig = make_subplots(rows=1, cols=3, 
                        specs=[[{'type':'pie'}, {'type':'pie'}, {'type':'pie'}]],
                        subplot_titles=("WT-IMP2", "mCherry-IMP2", "WT-mCherry"),
                        horizontal_spacing=0.1)
    fig.add_trace(px.pie(result_df_count[result_df_count['comparison_type'].str.contains('wt_imp2')],
                         labels='region', names='region', values='count').data[0], 
                  1, 1)
    fig.add_trace(px.pie(result_df_count[result_df_count['comparison_type'].str.contains('mcherry_imp2')],
                         labels='region', names='region', values='count').data[0], 
                  1, 2)
    fig.add_trace(px.pie(result_df_count[result_df_count['comparison_type'].str.contains('wt_mcherry')],
                         labels='region', names='region', values='count').data[0], 
                  1, 3)

    # Use `hole` to create a donut-like pie chart
    fig.update_traces(selector='pie', textinfo="value+percent",
                      texttemplate='%{value}<br>(%{percent})',
                      marker_colors=px.colors.qualitative.Pastel)
    fig.update_layout(legend=dict(orientation="h", yanchor="bottom", xanchor="center", x=0.5))

    fig.write_html(analyzed_result_path.joinpath('annotated_' + threshold + '.html'))
    fig.write_image(analyzed_result_path.joinpath('annotated_' + threshold + '.png'), scale=4)
    
    fig = px.sunburst(result_df_count, values='count', path=['comparison_type', 'region'], 
                      color='region', color_discrete_sequence=px.colors.qualitative.Pastel)
    fig.update_traces(textinfo='label+percent parent')
    fig.update_layout(uniformtext=dict(minsize=15, mode='show'))
    fig.write_html(analyzed_result_path.joinpath('annotated_combined_' + threshold + '.html'))
    fig.write_image(analyzed_result_path.joinpath('annotated_combined_'+ threshold + '.png'), scale=4)
    