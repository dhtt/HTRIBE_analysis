#conda env python3.11 
'''
This script gather all MEME-suite/XSTREME motives and annotate them to regions from a reference genome gtf file.

Usage: analyze_meme_motif.py 
'''

from pathlib import Path
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots

if __name__=='__main__':
    STREME_main_path = Path('/home/dhthutrang/HTRIBE_analysis/1%_bed_motif/XSTREME/_all/')
    motives_df = pd.read_csv(STREME_main_path.joinpath('motives_annotated.csv'), header=0, sep='\t')
    
#     motives_df_freq = motives_df \
#         .groupby(['comparison', 'motif', 'region']).size() \
#         .reset_index() 
#     motives_df_freq.columns = ['comparison_type', 'motif', 'region', 'count']
#     motives_df_freq['motif'] = motives_df_freq['motif'].apply(lambda x: x.split('-')[1])
#     motives_df_freq['perc'] = round(100*motives_df_freq['count']/ motives_df_freq.groupby(['comparison_type', 'motif'])['count'].transform('sum'), 2)
#     motives_df_freq.to_csv(STREME_main_path.joinpath('motives_annotated_freq.csv'), header=True, sep='\t')
    motives_df_freq = pd.read_csv(STREME_main_path.joinpath('motives_annotated_freq.csv'), header=0, sep='\t')
   
    fig = make_subplots(rows=1, cols=3, 
                        specs=[[{'type':'sunburst'}, {'type':'sunburst'}, {'type':'sunburst'}]],
                        subplot_titles=("WT-IMP2<br>vs.<br>WT-mCherry", "mCherry-IMP2<br>vs.<br>WT-mCherry", "WT-IMP2<br>vs.<br>mCherry-IMP2"),
                        horizontal_spacing=0.1)
    # fig.add_trace(px.pie(motives_df_freq[motives_df_freq['comparison_type'].str.contains('wt_imp2_vs_wt_mcherry')],
    #                     labels='region', names='region', values='count').data[0], 
    #                 1, 1)
    # fig.add_trace(px.pie(motives_df_freq[motives_df_freq['comparison_type'].str.contains('mcherry_imp2_vs_wt_mcherry')],
    #                     labels='region', names='region', values='count').data[0], 
    #                 1, 2)
    # fig.add_trace(px.pie(motives_df_freq[motives_df_freq['comparison_type'].str.contains('wt_imp2_vs_mcherry_imp2')],
    #                     labels='region', names='region', values='count').data[0], 
    #                 1, 3)

    # # Use `hole` to create a donut-like pie chart
    # fig.update_traces(selector='pie', textinfo="value+percent",
    #                     texttemplate='%{value}<br>(%{percent})',
    #                     marker_colors=px.colors.qualitative.Pastel)
    # fig.update_layout(legend=dict(orientation="h", yanchor="bottom", xanchor="center", x=0.5))

#     fig.write_html(STREME_main_path.joinpath('motives_annotated_pie.html'))
#     fig.write_image(STREME_main_path.joinpath('motives_annotated_pie.png'), scale=4)
    
    
    plot_colors = {'(?)':'#EDEDED', 'CDS':'#FFD966', '3UTR':'#97D077', '5UTR':'#EA6B66'}
    i = 1
    for part, col in plot_colors.items():
          fig.add_annotation(dict(font=dict(color=col, ), 
                                  text='<b>'+part+'</b>', 
                                  y=0, x=1/5*i,
                                  showarrow=False,
                                  xanchor='center'))
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
    
    # fig = px.sunburst(motives_df_freq, values='count', path=['comparison_type', 'motif', 'region'], 
    #                     color='region', color_discrete_sequence=px.colors.qualitative.Pastel)
    # fig.update_traces(textinfo='label+percent parent')
    # fig.update_layout(uniformtext=dict(minsize=15, mode='show'))
    
    
#     fig.update_traces(selector='sunburst', textinfo="percent parent")
    
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
    