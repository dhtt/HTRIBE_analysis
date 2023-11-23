from helper_functions import gather_read_counts, create_new_folder
from pathlib import Path
import pandas as pd
from sklearn.decomposition import PCA
import plotly.express as px
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os 

sns.set_theme(style="white")


if __name__=='__main__':
    # Define input and output paths
    quantified_transcripts_path = '/home/dhthutrang/TRIBE/mRNA_seq/processed/alignment_salmon_ncbi'
    analyzed_result_path = Path(quantified_transcripts_path).joinpath('analyzed_result')
    expression_data_path = Path(analyzed_result_path.joinpath('expression_data.csv'))
    create_new_folder(analyzed_result_path)
    
    # Get expression data
    # if os.path.exists(expression_data_path):
    #     expression = pd.read_csv(analyzed_result_path.joinpath('expression_data.csv'), index_col=0)
    # else:
    expression = gather_read_counts(quantified_transcripts_path, normalize=True)
    expression.to_csv(analyzed_result_path.joinpath('expression_data.csv'), sep='\t')
    
    # # PCA analysis 
    # pca = PCA(n_components=4)
    # pca_expression = pca.fit_transform(expression)
    # print(pca_expression)
    # pca_columns = ['PCA_' + str(i) for i in range(1, pca_expression.shape[1]+1, 1)]
    # pca_expression = pd.DataFrame(pca_expression, columns=pca_columns)
    # pca_expression['group'] = ['WT', 'WT', 'WT', 'mCherry', 'mCherry', 'mCherry', 'IMP2', 'IMP2', 'IMP2']
    # pca_expression['type'] = 'control'
    # pca_expression['type'][pca_expression['group'] == 'IMP2'] = 'test'   
    
    # ## Define plot params
    # var_ratio = pca.explained_variance_ratio_ 
    # total_var = var_ratio.sum() * 100
    # loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
    # labels = {'PCA_' + str(i+1): f"PC {i+1} ({var:.1f}%)" for i, var in enumerate(var_ratio * 100)}
    # marker_size = 12
    
    # ## Plot 3D PCA
    # fig = px.scatter_3d(pca_expression, x=pca_columns[0], y=pca_columns[1], z=pca_columns[2], 
    #                     color=pca_expression['group'], symbol=pca_expression['type'],
    #                     title=f'Total Explained Variance: {total_var:.2f}%', labels=labels,
    #                     color_discrete_sequence=px.colors.qualitative.Pastel)
    # fig.update_traces(marker_size=marker_size)
    # fig.write_html(analyzed_result_path.joinpath('PCA_3D.html'))
    
    # ## Plot 2D PCA
    # fig = px.scatter_matrix(pca_expression, labels=labels, dimensions=pca_columns, 
    #                         color=pca_expression['group'], symbol=pca_expression['type'],
    #                         title=f'Total Explained Variance: {total_var:.2f}%',
    #                         color_discrete_sequence=px.colors.qualitative.Pastel)
    # fig.update_traces(diagonal_visible=False, marker_size=marker_size)
    # fig.write_image(analyzed_result_path.joinpath('PCA_2D.png'), width=8, height=8)
    # fig.write_html(analyzed_result_path.joinpath('PCA_2D.html'))

    # Correlation analysis 
    sns.set(font_scale=1.2)
    correlation_types = ['pearson', 'spearman']
    for correlation_type in correlation_types:     
        count_df = expression.transpose()
        diag_mat = np.triu(np.ones_like(count_df.corr(correlation_type)))
        # cor_df_mat[diag_mat] = np.nan
        # cor_df = pd.DataFrame(cor_df_mat, columns=cor_df_colnames, index=cor_df_colnames)
        # print(cor_df)

        fig = sns.heatmap(count_df.corr(correlation_type), annot=True, cmap='rocket_r', 
                          mask=diag_mat, square=True, linewidths=.5, cbar_kws={"shrink": .5},
                          annot_kws = {'size': 12})
        plt.savefig(analyzed_result_path.joinpath('correlation_' + correlation_type + '.png'),
                    dpi=400)
        plt.close()
     