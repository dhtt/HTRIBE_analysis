#!/usr/bin/env python # [1]
"""
This script analyzes the HYPERTRIBE summarized results for each threshold generated (0%, 1%, 5%)

Usage: analyze_by_threshold.py HTRIBE_result_path threshold
"""

import sys
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
from matplotlib_venn import venn3_circles, venn3
from itertools import combinations


def read_bedgraph(result_dir: str, threshold: str):
    """This function reads HYPERTRIBE result files in xls format and return a data frame of concatenated results

    Args:
        result_dir (str): absolute path to the HYPERTRIBE result folder
        threshold (str): threshold of HYPERTRIBE analysis to be analyzed. This should be a string starting with _A2G.
            For example, "_A2G_1%"
    """
    all_files = os.listdir(result_dir)
    df_list = []
    for file in all_files:
        if file.endswith(threshold + ".xls"):
            comparison_type = "_".join(file.split('_')[1:3])
            df = pd.read_csv(os.path.join(result_dir, file), sep='\t', header=0, index_col=None)
            df["comparison_type"] = comparison_type
            df_list.append(df)
    df_list_concat = pd.concat(df_list)
    return(df_list_concat)


if __name__ == "__main__":
    # HTRIBE_result_path = sys.argv[1]
    # threshold = sys.argv[2]

    HTRIBE_result_path = "/home/dhthutrang/TRIBE/mRNA_seq/processed/extract.trim.align.dedup/test1/result/collapsed_replicates/comparisons/combined_result"
    threshold = "_A2G_1%"

    # Create a folder storing results
    analyzed_result_path = Path(HTRIBE_result_path).joinpath('analyzed_result')
    if os.path.exists(analyzed_result_path) == False:
        os.mkdir(analyzed_result_path)

    # # Gather a dictionary of resulting genes from the HYPERTRIBE analysis for a certain threshold
    # result_df = read_bedgraph(result_dir=HTRIBE_result_path, threshold=threshold)
    # genes_by_comparison_type = dict()

    # for comparison_type in set(result_df['comparison_type']):
    #     no_editing_sites_by_gene = dict()
    #     genes = set(result_df['Gene_name'][result_df['comparison_type'] == comparison_type])
    #     for gene in genes:
    #         no_editing_sites_by_gene[gene] = result_df['Num_edit_sites'][(result_df['Gene_name'] == gene) & 
    #                                                                      (result_df['comparison_type'] == comparison_type)].values[0]
    #     genes_by_comparison_type[comparison_type] = no_editing_sites_by_gene
    
    # with open(analyzed_result_path.joinpath('genes_by_comparison_type.pickle'), 'wb') as file:
    #     pickle.dump(genes_by_comparison_type, file, protocol=pickle.HIGHEST_PROTOCOL)

    # # Generate Venn diagram to compare the overlap between the detected genes for different comparisons
    # fig, ax = plt.subplots(figsize=(6, 6),  nrows=1, ncols=1)
    # venn3(subsets=[set(genes_by_comparison_type['wt_mcherry'].keys()), 
    #                set(genes_by_comparison_type['wt_imp2'].keys()), 
    #                set(genes_by_comparison_type['mcherry_imp2'].keys())],
    #       set_labels=['WT vs. mCherry', 'WT vs. IMP2', 'mCherry vs. IMP2'])
    # venn3_circles(subsets=[set(genes_by_comparison_type['wt_mcherry'].keys()), 
    #                        set(genes_by_comparison_type['wt_imp2'].keys()), 
    #                        set(genes_by_comparison_type['mcherry_imp2'].keys())],
    #               alpha=0.4, ax=ax)
    # fig.savefig(analyzed_result_path.joinpath("overlap" + threshold + ".png"))
    # plt.close(fig)
    
    
    with open(analyzed_result_path.joinpath('genes_by_comparison_type.pickle'), 'rb') as file:
        genes_by_comparison_type = pickle.load(file)
    comparison_type = ["wt_mcherry", "wt_imp2", "mcherry_imp2"]
    comparison_pair = list(combinations(comparison_type, 2))
    
    for pair in comparison_pair:
        comparison_1, comparison_2 = pair[0], pair[1]
        mutual_genes = set(list(genes_by_comparison_type[comparison_1].keys()) + list(genes_by_comparison_type[comparison_2].keys()))
        coords = dict()
        for gene in mutual_genes:
            if gene in genes_by_comparison_type[comparison_1].keys():
                x_coord = genes_by_comparison_type[comparison_1][gene]
            else: 
                x_coord = 0
                
            if gene in genes_by_comparison_type[comparison_2].keys():
                y_coord = genes_by_comparison_type[comparison_2][gene]
            else:
                y_coord = 0
            coords[gene] = (x_coord, y_coord)
        print(coords)