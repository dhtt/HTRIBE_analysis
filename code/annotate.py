from helper_functions import read_bedgraph, create_new_folder
from pathlib import Path
import sys
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
from itertools import combinations
import sys
import pandas as pd

# pio.templates.default = "plotly_white"

def annotate_by_threshold(HTRIBE_result_path, threshold, region="CDS_UTR"):
       '''
       Analyzes the HYPERTRIBE summarized results for a specific threshold and region.

       Parameters:
              HTRIBE_result_path (str): The path to the HYPERTRIBE result directory. Default: "mRNA_seq/processed/extract.trim.align.dedup/test1/result/diff_span/all_span"
              threshold (str): The threshold value. Default: "1%"
              region (str, optional): The region to analyze. Default: "CDS_UTR".

       Returns:
              None
       '''
       workdir = HTRIBE_result_path

       # Create a folder storing results
       HTRIBE_result_path = workdir + '/' + region + '/OR/OR'
       analyzed_result_path = Path(HTRIBE_result_path).joinpath('analyzed_result')
       create_new_folder(analyzed_result_path)

       result_df = read_bedgraph(result_dir=HTRIBE_result_path, threshold=threshold, file_extension='.bedgraph', header=None)
       result_df.columns = ['chr', 'start', 'end', 'region', 'strand', 'number', 'comparison_type']
       # result_df['region'] = result_df['region'].str.split(pat=',')
       # Uncomment here to remove CDS if a region is also UTR
       # result_df['region'] = result_df['region'].apply(lambda x: list(filter(lambda i: i!="CDS", x)) if len(x) > 1 else x)
       # result_df = result_df.explode('region')
       result_df_count = result_df.groupby(['comparison_type', 'region'], as_index=False).count().iloc[:, :3]
       result_df_count.columns = ['comparison_type', 'region', 'count']
       result_df_count['perc'] = round(100*result_df_count['count']/sum(result_df_count['count']), 2)
       result_df_count['group_perc'] = round(100*result_df_count['count']/result_df_count.groupby('comparison_type', as_index=False)['count'].transform('sum')['count'], 2)
       result_df_count['comparison_type'] = pd.Categorical(
              result_df_count['comparison_type'],
              categories=['wt_imp2', 'mcherry_imp2', 'wt_mcherry', 'subwt_submcherry']
              )
       result_df_count['region'] = ['/'.join(x.split(',')) for x in result_df_count['region']]
       result_df_count = result_df_count.sort_values(['comparison_type', 'region'])
       
       # Make a bar plot
       fig = px.bar(result_df_count, color="region", y="group_perc", x="comparison_type",
              color_discrete_sequence=px.colors.qualitative.Set3, 
              text=['{}<br>({}%)'.format(v, p) for v, p in zip(result_df_count['count'], result_df_count['group_perc'])],
              labels={'group_perc': '<b>Percentage (%)</b>', 'comparison_type': '', 'region': '<b>Region</b>'}
              )
       fig.update_layout(
              uniformtext_minsize=13, uniformtext_mode='hide', font_family="Arial", font_color='black',
              xaxis = dict(tickangle=-45, tickmode = 'array', 
                     tickvals = ['wt_imp2', 'mcherry_imp2', 'wt_mcherry', 'subwt_submcherry'],
                     ticktext = ['<b>' + s + '</b>' for s in ['WT vs IMP2', 'mCherry vs IMP2', 'WT vs mCherry', 'WT/mCherry vs IMP2<br>-<br>WT vs mCherry']]
                     )
              )
       fig.update_traces(insidetextanchor="middle")
       fig.write_html(analyzed_result_path.joinpath('annotated_bar_' + threshold + '.html'))
       fig.write_image(analyzed_result_path.joinpath('annotated_bar_' + threshold + '.png'), scale=5, width=600, height=620)
       

       # Make a pie chart
       fig = make_subplots(rows=2, cols=2, 
                                          specs=[[{'type':'pie'}, {'type':'pie'}], [{'type':'pie'}, {'type':'pie'}]],
                                          subplot_titles=("WT-IMP2", "mCherry-IMP2","", "WT-mCherry", "WT/mCherry-IMP2"),
                                          horizontal_spacing=0.1)
       fig.add_trace(px.pie(result_df_count[result_df_count['comparison_type'].str.contains('wt_imp2')],
                                           labels='region', names='region', values='count').data[0], 
                              1, 1)
       fig.add_trace(px.pie(result_df_count[result_df_count['comparison_type'].str.contains('mcherry_imp2')],
                                           labels='region', names='region', values='count').data[0], 
                              1, 2)
       fig.add_trace(px.pie(result_df_count[result_df_count['comparison_type'].str.contains('wt_mcherry')],
                                           labels='region', names='region', values='count').data[0], 
                              2, 1)
       fig.add_trace(px.pie(result_df_count[result_df_count['comparison_type'].str.contains('subwt_submcherry')],
                                           labels='region', names='region', values='count').data[0], 
                              2, 2)
       

       # Make a donut-like pie chart
       fig.update_traces(selector='pie', textinfo="value+percent",
                                     texttemplate='%{value}<br>(%{percent})',
                                     marker_colors=px.colors.qualitative.Pastel)
       fig.update_layout(legend=dict(orientation="h", yanchor="bottom", xanchor="center", x=0.5))

       print("Saving the annotated results to: " + str(analyzed_result_path))

       fig.write_html(analyzed_result_path.joinpath('annotated_' + threshold + '.html'))
       fig.write_image(analyzed_result_path.joinpath('annotated_' + threshold + '.png'), scale=4, width=1200, height=1000)
       
       fig = px.sunburst(result_df_count, values='count', path=['comparison_type', 'region'], 
                                     color='region', color_discrete_sequence=px.colors.qualitative.Pastel)
       fig.update_traces(textinfo='label+percent parent')
       fig.update_layout(uniformtext=dict(minsize=15, mode='show'))
       fig.write_html(analyzed_result_path.joinpath('annotated_combined_' + threshold + '.html'))
       fig.write_image(analyzed_result_path.joinpath('annotated_combined_'+ threshold + '.png'), scale=4, width=1200, height=1000)

if __name__=='__main__':
       annotate_by_threshold(sys.argv[1], sys.argv[2], sys.argv[3])
    