#!/usr/bin/env python # [1]
'''
This script analyzes the HYPERTRIBE summarized results for each threshold generated (0%, 1%, 5%)

Usage: analyze_by_threshold.py HTRIBE_result_path threshold
'''

import seaborn as sns
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
from matplotlib_venn import venn3_circles, venn3
from itertools import combinations
from math import nan


def read_bedgraph(result_dir: str, threshold: str):
    '''This function reads HYPERTRIBE result files in xls format and return a data frame of concatenated results

    Args:
        result_dir (str): absolute path to the HYPERTRIBE result folder
        threshold (str): threshold of HYPERTRIBE analysis to be analyzed. This should be a string starting with _A2G.
            For example, '_A2G_1%'
    '''
    all_files = os.listdir(result_dir)
    df_list = []
    for file in all_files:
        if file.endswith(threshold + '.xls'):
            comparison_type = '_'.join(file.split('_')[1:3])
            df = pd.read_csv(os.path.join(result_dir, file),
                             sep='\t', header=0, index_col=None)
            df['comparison_type'] = comparison_type
            df_list.append(df)
    df_list_concat = pd.concat(df_list)
    return(df_list_concat)


def get_sorted_labels(ax_):
    handlers, labels = ax.get_legend_handles_labels()
    handlers_dict = dict(zip(labels, handlers))
    handlers_dict = dict(sorted(handlers_dict.items()))
    return handlers_dict


def generate_venn_diagram(genes_by_comparison_type_, analyzed_result_path_, threshold_):
    fig, ax = plt.subplots(figsize=(6, 6),  nrows=1, ncols=1)
    venn3(subsets=[set(genes_by_comparison_type_['wt_mcherry'].keys()),
                   set(genes_by_comparison_type_['wt_imp2'].keys()),
                   set(genes_by_comparison_type_['mcherry_imp2'].keys())],
          set_labels=['WT vs. mCherry', 'WT vs. IMP2', 'mCherry vs. IMP2'])
    venn3_circles(subsets=[set(genes_by_comparison_type_['wt_mcherry'].keys()),
                           set(genes_by_comparison_type_['wt_imp2'].keys()),
                           set(genes_by_comparison_type_['mcherry_imp2'].keys())],
                  alpha=0.4, ax=ax)
    fig.savefig(analyzed_result_path_.joinpath(
        'overlap_' + threshold_ + '.png'))
    plt.close(fig)


if __name__ == '__main__':
    # HTRIBE_result_path = sys.argv[1]
    # threshold = sys.argv[2]

    HTRIBE_result_path = '/home/dhthutrang/TRIBE/mRNA_seq/processed/extract.trim.align.dedup/test1/result/collapsed_replicates/comparisons/combined_result'
    threshold = '5%'
    comparison_type = ['wt_mcherry', 'wt_imp2', 'mcherry_imp2']
    comparison_pair = list(combinations(comparison_type, 2))

    # Create a folder storing results
    analyzed_result_path = Path(HTRIBE_result_path).joinpath('analyzed_result')
    if os.path.exists(analyzed_result_path) == False:
        os.mkdir(analyzed_result_path)

    # Gather a dictionary of resulting genes from the HYPERTRIBE analysis for a certain threshold.
    # Return genes_by_comparison_type[comparison_type][gene]: Number of editing sites
    result_df = read_bedgraph(
        result_dir=HTRIBE_result_path, threshold=threshold)
    genes_by_comparison_type = dict()

    for comparison_type in set(result_df['comparison_type']):
        no_editing_sites_by_gene = dict()
        genes = set(result_df['Gene_name']
                    [result_df['comparison_type'] == comparison_type])
        for gene in genes:
            no_editing_sites_by_gene[gene] = result_df['Num_edit_sites'][(result_df['Gene_name'] == gene) &
                                                                         (result_df['comparison_type'] == comparison_type)].values[0]
        genes_by_comparison_type[comparison_type] = no_editing_sites_by_gene

    with open(analyzed_result_path.joinpath('genes_by_comparison_type_' + threshold + '.pickle'), 'wb') as file:
        pickle.dump(genes_by_comparison_type, file,
                    protocol=pickle.HIGHEST_PROTOCOL)

    # Generate Venn diagram to compare the overlap between the detected genes for different comparisons
    generate_venn_diagram(genes_by_comparison_type_=genes_by_comparison_type,
                          analyzed_result_path_=analyzed_result_path, threshold_=threshold)

    # Store number of editing sites identified for each gene from a pairwise comparison in coords_by_no_editing_sites
    # Return coords_by_no_editing_sites[comparison_pair][gene]: (x, y) while x is the number of editing sites for gene in comparison 1
    # and y is the number of editing site for gene in comparison 2)
    with open(analyzed_result_path.joinpath('genes_by_comparison_type_' + threshold + '.pickle'), 'rb') as file:
        genes_by_comparison_type = pickle.load(file)
    coords_by_no_editing_sites = dict()

    # For each comparison (wild type vs. mCherry or IMP2 or mCherry vs. IMP2)
    for pair in comparison_pair:
        comparison_1, comparison_2 = pair[0], pair[1]
        mutual_genes = set(list(genes_by_comparison_type[comparison_1].keys(
        )) + list(genes_by_comparison_type[comparison_2].keys()))
        coords = dict()

        # Store the number of editing sites from comparison 1 as x-coordinate and from 2 as y-coordinate
        for gene in mutual_genes:
            if gene in genes_by_comparison_type[comparison_1].keys():
                x_coord = genes_by_comparison_type[comparison_1][gene]
            else:
                x_coord = nan

            if gene in genes_by_comparison_type[comparison_2].keys():
                y_coord = genes_by_comparison_type[comparison_2][gene]
            else:
                y_coord = nan
            coords[gene] = (x_coord, y_coord)
        coords_by_no_editing_sites[str(
            comparison_1 + ':' + comparison_2)] = coords

    with open(analyzed_result_path.joinpath('coords_by_no_editing_sites' + threshold + '.pickle'), 'wb') as file:
        pickle.dump(coords_by_no_editing_sites, file,
                    protocol=pickle.HIGHEST_PROTOCOL)

    # Create a dataframe for scatter plot comparing the gene result identified from the comparisons
    with open(analyzed_result_path.joinpath('coords_by_no_editing_sites' + threshold + '.pickle'), 'rb') as file:
        coords_by_no_editing_sites = pickle.load(file)
    comparison_pair = list(coords_by_no_editing_sites.keys())
    all_comparisons = []
    for pair in comparison_pair:
        comparison_1 = pair.split(':')
        comparison_2 = comparison_1[1]
        comparison_1 = comparison_1[0]

        df = pd.DataFrame.from_dict(
            coords_by_no_editing_sites[pair]).transpose()
        df.columns = ['comparison_1', 'comparison_2']
        df['presence'] = 'Both comparisons'
        df['presence'][df[['comparison_1', 'comparison_2']].isnull().any(1)] = 'Either comparison'
        df = df.fillna(0)
        df['comparison_type'] = comparison_1 + '_' + comparison_2
        df['gene_name'] = df.index
        axis_lim = max(max(df['comparison_1']), max(df['comparison_2']))
        
        # Plot genes identified from the comparisons by the number of editing sites
        fig, ax = plt.subplots(figsize=(5, 4),  nrows=1, ncols=1)
        ax = sns.scatterplot(data=df, x='comparison_1', y='comparison_2',
                             hue='presence', palette={'Both comparisons': 'green', 'Either comparison': 'orange'},
                             style='presence', markers={'Both comparisons': 'o', 'Either comparison': 'X'},
                             linewidth=0, alpha=0.4)
        ax.plot([0, axis_lim], [0, axis_lim], 'red', linestyle='dashed', alpha = 0.6)

        ax.set(title='Genes with A2G identified from the\npairwise comparisons (T = ' + threshold + ')',
               xlabel=comparison_1, ylabel=comparison_2)

        handlers_dict = get_sorted_labels(ax)
        plt.legend(handlers_dict.values(), handlers_dict.keys(),
                   title='Genes with A2G sites\nidentified in:')
        fig.savefig(analyzed_result_path.joinpath(
            'scatterplot_' + comparison_1 + ':' + comparison_2 + '_' + threshold + '.png'))
        plt.close(fig)

        # Store the dataframe for this comparison pair to make a combined plot
        all_comparisons.append(df)

    # Combined scatter plot from all type of pair comparisons at once
    all_comparisons_df = pd.concat(all_comparisons)
    special_genes = all_comparisons_df[(all_comparisons_df['presence'] == 'Both comparisons') & (all_comparisons_df['comparison_type'] == 'wt_imp2_mcherry_imp2')]
    special_genes = special_genes.sort_values(by=['comparison_1', 'comparison_2'])
    print(special_genes.tail(10))
    
    axis_lim = max(max(all_comparisons_df['comparison_1']), max(all_comparisons_df['comparison_2']))
    fig, ax = plt.subplots(figsize=(5, 4),  nrows=1, ncols=1)
    ax = sns.scatterplot(data=all_comparisons_df, x='comparison_1', y='comparison_2',
                         hue='comparison_type', alpha=0.4, linewidth=0.5)
    ax.plot([0, axis_lim], [0, axis_lim], 'red', linestyle='dashed', alpha = 0.6)
    ax.set(title='Genes with A2G identified from the\npairwise comparisons (T = ' + threshold + ')',
           xlabel='Comparison 1', ylabel='Comparison 2')
    plt.legend(title='Comparison type',
               labels=['_', 'WT_mCherry vs. WT_IMP2', 'WT_mCherry vs. mCherry_IMP2', 'WT_IMP2 vs. WT_mCherry'])
    fig.savefig(analyzed_result_path.joinpath(
        'scatterplot_' + threshold + '.png'))
    plt.close(fig)
