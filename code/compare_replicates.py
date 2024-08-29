"""
This script analyzes the replicates and computes the Jaccard similarity coefficient for filtered pairs.

Usage: python compare_replicates.py <replicate_result_path>

Command line arguments:
- replicate_result_path (str): The path to the replicate result directory. Default: "~/TRIBE/mRNA_seq/processed/extract.trim.align.dedup/test1/result/compare_replicates"

"""

import os
import sys
import pickle
from pathlib import Path
from itertools import combinations
from helper_functions import read_bedgraph

def compute_jaccard(set_1, set_2):
    """
    Computes the Jaccard similarity coefficient between two sets.

    Parameters:
    set_1 (set): The first set.
    set_2 (set): The second set.

    Returns:
    float: The Jaccard similarity coefficient.
    """
    no_intersection = len(set(set_1).intersection(set(set_2)))
    no_union = len(set(set_1).union(set(set_2)))
    return(no_intersection/no_union)

def analyze_replicates(replicate_result_path, threshold):
    """
    Analyzes the replicates and computes the Jaccard similarity coefficient for filtered pairs.

    Parameters:
    replicate_result_path (str): The path to the replicate result directory.
    threshold (str): The threshold for the analysis.

    Returns:
    None
    """
    analyzed_result_path = Path(replicate_result_path).joinpath('analyzed_result')

    if os.path.exists(analyzed_result_path) == False:
        os.mkdir(analyzed_result_path)
    
    genes_by_comparison_type_filename = analyzed_result_path.joinpath('genes_by_comparison_type_' + threshold + '.pickle')
    if os.path.exists(genes_by_comparison_type_filename):
        with open(genes_by_comparison_type_filename, 'rb') as file:
            genes_by_comparison_type = pickle.load(file)
    else:
        # Gather a dictionary of resulting genes from the HYPERTRIBE analysis for a certain threshold.
        # Return genes_by_comparison_type[comparison_type][gene]: Number of editing sites
        result_df = read_bedgraph(result_dir=replicate_result_path, threshold=threshold)
        genes_by_comparison_type = dict()

        for comparison_type in set(result_df['comparison_type']):
            no_editing_sites_by_gene = dict()
            genes = set(result_df['Gene_name'][result_df['comparison_type'] == comparison_type])
            for gene in genes:
                no_editing_sites_by_gene[gene] = result_df['Num_edit_sites'][(result_df['Gene_name'] == gene) &
                                                                            (result_df['comparison_type'] == comparison_type)].values[0]
            genes_by_comparison_type[comparison_type] = no_editing_sites_by_gene
        
        with open(analyzed_result_path.joinpath('genes_by_comparison_type_' + threshold + '.pickle'), 'wb') as file:
            pickle.dump(genes_by_comparison_type, file,
                        protocol=pickle.HIGHEST_PROTOCOL)
            
    
    pairs = list(genes_by_comparison_type.keys())
    pairs.sort()
    compared_pairs = list(combinations(pairs, 2))
    compared_pairs = sorted(compared_pairs, key = lambda elem: (elem[0], elem[1]))
    
    filtered_pairs = [pair for pair in compared_pairs if ("_7" in pair[0] and "_7" in pair[1]) or ("_8" in pair[0] and "_8" in pair[1]) or ("_9" in pair[0] and "_9" in pair[1])]

    for pair in filtered_pairs:
        jaccard_index = compute_jaccard(genes_by_comparison_type[pair[0]].keys(), genes_by_comparison_type[pair[1]].keys())
        print(f"{pair}: {round(jaccard_index, 2)}")

if __name__ == '__main__':
    replicate_result_path = sys.argv[1]
    threshold = '5%'
    wt = [1, 2, 3]
    mcherry = [4, 5, 6]
    imp2 = [7, 8, 9]

    # Create a folder storing results
    analyze_replicates(replicate_result_pathreplicate_result_path=replicate_result_path, threshold=threshold)
