"""
This script performs PCA analysis and correlation analysis on gene expression data.

The script reads expression data from a CSV file or computes it from raw read counts.
It then performs PCA analysis on the expression data and generates 2D and 3D scatter plots.
The total explained variance is calculated and displayed in the plots.

The script also performs correlation analysis on the expression data and generates a heatmap.

The generated plots and heatmap are saved as image files.

"""

from helper_functions import gather_read_counts, create_new_folder
from pathlib import Path
import pandas as pd
from sklearn.decomposition import PCA
import plotly.express as px
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import os 


if __name__=='__main__':
    # Define input and output paths
    quantified_transcripts_path = '/home/dhthutrang/TRIBE/mRNA_seq/processed/alignment_salmon_ncbi'
    analyzed_result_path = Path(quantified_transcripts_path).joinpath('analyzed_result')
    expression_data_path = Path(analyzed_result_path.joinpath('expression_data.csv'))
    create_new_folder(analyzed_result_path)
    
    # Get expression data
    if os.path.exists(expression_data_path):
        expression = pd.read_csv(analyzed_result_path.joinpath('expression_data.csv'), index_col=0, sep='\t')
        expression.index = [" ".join(x.split("_")) for x in expression.index]
    else:
        expression = gather_read_counts(quantified_transcripts_path, normalize=True)
        expression.to_csv(analyzed_result_path.joinpath('expression_data.csv'), sep='\t')
    
    # PCA analysis 
    pca = PCA(n_components=4)
    pca_expression = pca.fit_transform(expression)
    pca_columns = ['PCA_' + str(i) for i in range(1, pca_expression.shape[1]+1, 1)]
    pca_expression = pd.DataFrame(pca_expression, columns=pca_columns)
    pca_expression['group'] = ['WT', 'WT', 'WT', 'mCherry', 'mCherry', 'mCherry', 'IMP2', 'IMP2', 'IMP2']
    pca_expression['type'] = 'Control'
    pca_expression['type'][pca_expression['group'] == 'IMP2'] = 'Test'   
    
    # Define plot params
    var_ratio = pca.explained_variance_ratio_ 
    total_var = var_ratio.sum() * 100
    loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    labels = {'PCA_' + str(i+1): f"PC {i+1} ({var:.1f}%)" for i, var in enumerate(var_ratio * 100)}
    marker_size = 15
    
    # Plot 3D PCA
    fig = px.scatter_3d(pca_expression, x=pca_columns[0], y=pca_columns[1], z=pca_columns[2], 
                        color=pca_expression['group'], symbol=pca_expression['type'], 
                        title=f'Total Explained Variance: {total_var:.2f}%', labels=labels,
                        color_discrete_map={'WT': '#97D077', 'mCherry': '#A680B8', 'IMP2': '#FFD966'})
    fig.update_layout(template='seaborn', margin=dict(l=0, r=0, b=0, t=0), 
                      scene_aspectmode='manual', scene_aspectratio=dict(x=1.5, y=1, z=1),
                      legend=dict(orientation="h", xanchor='center', x=0.5),
                      legend_title_text='Group / Sample type')
    fig.update_traces(marker=dict(size=marker_size, line=dict(color='black', width=2)),
                      selector=dict(mode='markers'))
    fig.write_image(analyzed_result_path.joinpath('PCA_2D.png'), width=800, height=800, scale=4)
    fig.write_html(analyzed_result_path.joinpath('PCA_3D.html'))
    
    ## Plot 2D PCA
    fig = px.scatter_matrix(pca_expression, labels=labels, dimensions=pca_columns, 
                            color=pca_expression['group'], symbol=pca_expression['type'],
                            title=f'Principle components analysis of expression data (Total Explained Variance: {total_var:.2f}%)', opacity=0.9,
                            color_discrete_map={'WT': '#97D077', 'mCherry': '#A680B8', 'IMP2': '#FFD966'})
    fig.update_layout(template='seaborn', margin_pad=2, 
                      legend=dict(orientation="h", xanchor='center', x=0.5, yanchor='bottom', y=1.02),
                      legend_title_text='Group / Sample type',
                      margin=dict(l=100)
                      )
    fig.update_traces(diagonal_visible=False, marker=dict(size=marker_size-3, line=dict(color='black', width=2)))
    fig.write_image(analyzed_result_path.joinpath('PCA_2D.png'), width=800, height=600, scale=4)
    fig.write_html(analyzed_result_path.joinpath('PCA_2D.html'))

    # Correlation analysis 
    sns.set_theme(style="white")
    cmap = cm.get_cmap('spring_r') # Get your favorite cmap
    new_cmap = ListedColormap(cmap(np.linspace(0, 0.25, 256)))
    
    correlation_data = {"pearson": tuple(), "spearman": tuple()}
    for correlation_type in correlation_data.keys():     
        count_df = expression.transpose()
        if correlation_type == "pearson":
            diag_mat = np.triu(np.ones_like(count_df.corr(correlation_type)))
        elif correlation_type == "spearman":
            diag_mat = np.tril(np.ones_like(count_df.corr(correlation_type)))
        correlation_data[correlation_type] = (count_df, diag_mat)
        
    fig = plt.gcf()
    fig.set_size_inches(8, 8)
    fig = sns.heatmap(correlation_data['pearson'][0].corr('pearson'), mask=correlation_data['pearson'][1], 
                      annot=True, cmap='summer_r', square=True, linewidths=.5, 
                      cbar_kws={'label': 'Pearson', 'location': 'right', "shrink": .5}, annot_kws = {'size': 12})
    fig = sns.heatmap(correlation_data['spearman'][0].corr('spearman'), mask=correlation_data['spearman'][1], 
                      annot=True, cmap=new_cmap, square=True, linewidths=.5, 
                      cbar_kws={'label': 'Spearman', 'location': 'right', "shrink": .5}, annot_kws = {'size': 12})
    plt.title("Gene expression correlation", fontdict={'weight': 'bold', 'size': 15}, pad=15)
    plt.savefig(analyzed_result_path.joinpath('correlation_' + 'pearson_spearman' + '.png'),
                dpi=400,  bbox_inches='tight')
    plt.close()