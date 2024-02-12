'''
This script analyzes the HYPERTRIBE summarized results for each threshold generated (0%, 1%, 5%)

Usage: analyze_by_threshold.py HTRIBE_result_path threshold region

This script performs the following tasks:
1. Reads the HYPERTRIBE summarized results from a specified directory.
2. Generates Venn diagrams to compare the overlap between the detected genes for different comparisons.
3. Reports significant genes and their comparison pairs.
4. Computes the Jaccard distance based on shared genes between comparison pairs.
5. Compares the analysis results to Jeet's results.
6. Plots common genes between comparison pairs in a scatter plot.

The script takes three command-line arguments:
- HTRIBE_result_path: The path to the HYPERTRIBE result directory.
- threshold: The threshold value for the analysis.
- region: If region is blank, the script will analyze all regions. Otherwise, it will analyze the specified region.

'''

from helper_functions import read_bedgraph
import seaborn as sns
import os
import sys
import pickle
import pandas as pd
import string
import matplotlib.pyplot as plt
from pathlib import Path
from collections import defaultdict
from matplotlib_venn import venn3_circles, venn3
from itertools import combinations
from math import nan
from helper_functions import create_new_folder
from compare_replicates import compute_jaccard
plt.style.use('seaborn-v0_8-bright')
plt.rc('font', family='Arial') 

def get_sorted_labels(ax_):
    """
    Returns a dictionary of legend labels sorted alphabetically.

    Parameters:
    - ax_ (matplotlib.axes.Axes): The axes object containing the legend.

    Returns:
    - dict: A dictionary mapping legend labels to their corresponding handlers, sorted alphabetically.
    """
    handlers, labels = ax_.get_legend_handles_labels()
    handlers_dict = dict(zip(labels, handlers))
    handlers_dict = dict(sorted(handlers_dict.items()))
    return handlers_dict


def generate_venn_diagram(genes_by_comparison_type_, analyzed_result_path_, threshold_):
    """
    Generates a Venn diagram based on the given genes by comparison type.

    Parameters:
    - genes_by_comparison_type_ (dict): A dictionary containing genes categorized by comparison type.
    - analyzed_result_path_ (str): The path where the analyzed result will be saved.
    - threshold_ (str): The threshold value used for analysis.

    Returns:
    - None
    """
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


def gather_genes_by_comparison_type(analyzed_result_path, HTRIBE_result_path, threshold):
    """
    Gathers genes by comparison type from the HYPERTRIBE analysis results for a certain threshold.

    Parameters:
    - analyzed_result_path (str): The path where the analyzed result will be saved.
    - HTRIBE_result_path (str): The path to the HYPERTRIBE result directory.
    - threshold (str): The threshold value used for analysis.

    Returns:
    - genes_by_comparison_type (dict): A dictionary containing genes categorized by comparison type.
    """

    if os.path.exists(analyzed_result_path.joinpath('genes_by_comparison_type_' + threshold + '.pickle')):
        with open(analyzed_result_path.joinpath('genes_by_comparison_type_' + threshold + '.pickle'), 'rb') as file:
            genes_by_comparison_type = pickle.load(file)
    else:
        result_df = read_bedgraph(result_dir=HTRIBE_result_path, threshold=threshold, file_extension='.xls')
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

    return genes_by_comparison_type
   
   
def report_sig_genes(sig_gene_filename, threshold, genes_by_comparison_type):
    """
    Writes the significant genes to a file.

    Args:
        - sig_gene_filename (str): The path to the output file.
        - threshold (str): The threshold value.
        - genes_by_comparison_type (dict): A dictionary containing the significant genes for each comparison pair.

    Returns:
        - None
    """
    if os.path.exists(sig_gene_filename):
        pass
    else:
        with open(sig_gene_filename, 'a') as outfile:
            outfile.write('Threshold: ' + threshold + '\n')

        with open(sig_gene_filename, 'a') as outfile:
            for comparison_pair, sig_genes in genes_by_comparison_type.items():
                outfile.write("*" * 30 + '\n')
                outfile.write(comparison_pair + '\n')
                outfile.write('","'.join(list(sig_genes.keys())))
                

def write_pairwise_jaccard_distance(main_output_filename, threshold, genes_by_comparison_type, comparison_pair):
    """
    Writes the pairwise Jaccard distance between genes to a file.

    Parameters:
    - main_output_filename (str): The filename of the output file.
    - threshold (str): The threshold value.
    - genes_by_comparison_type (dict): A dictionary containing genes grouped by comparison type.
    - comparison_pair (list): A list of comparison pairs.

    Returns:
    - None
    """
    
    if os.path.exists(main_output_filename):
        pass
    else:
        with open(main_output_filename, 'a') as outfile:
            outfile.write('## PAIRWISE JACCARD DISTANCE ##\nThreshold: ' + threshold + '\n')
            
        for pair in comparison_pair:
            result = pair[0] + '_' + pair[1] + ": " + str(round(compute_jaccard(list(genes_by_comparison_type[pair[0]].keys()), list(genes_by_comparison_type[pair[1]].keys())), 2))
            print(result)
            with open(main_output_filename, 'a') as outfile:
                outfile.write(result + '\n')

def compare_with_jeet_result(Jeet_result_path, threshold, genes_by_comparison_type, comparison_type, main_output_filename):
    """
    Compare analysis results to Jeet's results.

    Parameters:
    - Jeet_result_path (str): The path to Jeet's result directory.
    - threshold (float): The threshold value.
    - genes_by_comparison_type (dict): A dictionary containing genes grouped by comparison type.
    - comparison_type (list): A list of comparison types.
    - main_output_filename (str): The filename of the main output file.

    Returns:
    - None
    """
    # Compare analysis results to Jeet's results
    Jeet_result = read_bedgraph(result_dir=Jeet_result_path, threshold=threshold, file_extension='_results.xls')
    Jeet_result_genes = list(Jeet_result['Gene_name'])
    gene_list = [gene for gene in list(genes_by_comparison_type['mcherry_imp2'].keys())]
    
    for comp_type in comparison_type:
        gene_list = list(genes_by_comparison_type[comp_type].keys())
        no_intersecting_genes = len(set(gene_list).intersection(set(Jeet_result_genes)))
        result = str(round(compute_jaccard(set(gene_list), set(Jeet_result_genes))*1000, 2))
        with open(main_output_filename, 'a') as outfile:
            outfile.write("## JACCARD DISTANCE TO JEET'S RESULT ##\n" +
                        comp_type + ': ' + result + '-' + str(no_intersecting_genes) + '\n')
            

def generate_coords_by_no_editing_sites(coords_by_no_editing_sites_filename, comparison_pair, genes_by_comparison_type):
    """
    Generate coordinates by no editing sites.

    Parameters:
    coords_by_no_editing_sites_filename (str): The filename to save the coords_by_no_editing_sites dictionary.
    comparison_pair (list): The list of comparison pairs.
    genes_by_comparison_type (dict): The dictionary of genes by comparison type.

    Returns:
    dict: The coords_by_no_editing_sites dictionary.
    """
    
    if os.path.exists(coords_by_no_editing_sites_filename):
        # Load the existing coords_by_no_editing_sites dictionary from file
        with open(coords_by_no_editing_sites_filename, 'rb') as file:
            coords_by_no_editing_sites = pickle.load(file)
    else:
        # Create a new coords_by_no_editing_sites dictionary
        coords_by_no_editing_sites = dict()

        # For each comparison pair
        for pair in comparison_pair:
            comparison_1, comparison_2 = pair[0], pair[1]
            mutual_genes = set(list(genes_by_comparison_type[comparison_1].keys()) +
                                list(genes_by_comparison_type[comparison_2].keys()))
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
            coords_by_no_editing_sites[str(comparison_1 + ':' + comparison_2)] = coords

        # Save the coords_by_no_editing_sites dictionary to file
        with open(coords_by_no_editing_sites_filename, 'wb') as file:
            pickle.dump(coords_by_no_editing_sites, file, protocol=pickle.HIGHEST_PROTOCOL)

    return coords_by_no_editing_sites


def plot_scatterplots(coords_by_no_editing_sites, analyzed_result_path, threshold):
    """
    Plots scatterplots for pairwise comparisons of gene editing sites.

    Parameters:
    - coords_by_no_editing_sites (dict): Dictionary containing coordinates of gene editing sites for each pairwise comparison.
    - analyzed_result_path (str): Path to save the analyzed results.
    - threshold (str): Threshold value for the scatterplots.

    Returns:
    - None
    """
    
    all_comparisons = []
    comparison_pair = list(coords_by_no_editing_sites.keys())
    for i, pair in enumerate(comparison_pair):
        comparison_1 = pair.upper().replace("MCHERRY", "mCherry").split(':')
        comparison_2 = " vs. ".join(comparison_1[1].split("_")) 
        comparison_1 = " vs. ".join(comparison_1[0].split("_"))  

        df = pd.DataFrame.from_dict(coords_by_no_editing_sites[pair]).transpose()
        df.columns = ['comparison_1', 'comparison_2']
        df['presence'] = 'Both comparisons'
        df['presence'][df[['comparison_1', 'comparison_2']].isnull().any(1)] = 'Either comparison'
        df = df.fillna(0)
        df['comparison_type'] = comparison_1 + '_' + comparison_2
        df['gene_name'] = df.index
        axis_lim = max(max(df['comparison_1']), max(df['comparison_2']))
        
        # Plot genes identified from the comparisons by the number of editing sites
        fig, ax = plt.subplots(figsize=(5, 4), nrows=1, ncols=1)
        ax = sns.scatterplot(data=df, x='comparison_1', y='comparison_2',
                                hue='presence', style='presence', 
                                markers={'Both comparisons': 'o', 'Either comparison': 'X'},
                                palette={'Both comparisons': '#7dc754', 'Either comparison': '#faa719'},
                                linewidth=0.5, alpha=0.4)
        ax.plot([0, axis_lim], [0, axis_lim], 'red', linestyle='dashed', alpha = 0.6)
        ax.set(xlabel=comparison_1, ylabel=comparison_2)

        handlers_dict = get_sorted_labels(ax)
        plt.legend(handlers_dict.values(), handlers_dict.keys())
        plt.title(string.ascii_uppercase[i] + '. ' + comparison_1 + " against " + comparison_2,
                    fontweight='bold', loc='left')
        plt.gcf().set_size_inches(4, 4)
        fig.savefig(analyzed_result_path.joinpath('scatterplot_' + comparison_1 + ':' + comparison_2 + '_' + threshold + '.png'),
                    dpi=400, bbox_inches='tight')
        plt.close(fig)

        # Store the dataframe for this comparison pair to make a combined plot
        all_comparisons.append(df)

    # Combined scatter plot from all type of pair comparisons at once
    all_comparisons_df = pd.concat(all_comparisons)
    special_genes = all_comparisons_df[(all_comparisons_df['presence'] == 'Both comparisons') & (all_comparisons_df['comparison_type'] == 'wt_imp2_mcherry_imp2')]
    special_genes = special_genes.sort_values(by=['comparison_1', 'comparison_2'])
    
    axis_lim = max(max(all_comparisons_df['comparison_1']), max(all_comparisons_df['comparison_2']))
    fig, ax = plt.subplots(figsize=(5, 4),  nrows=1, ncols=1)
    ax = sns.scatterplot(data=all_comparisons_df, x='comparison_1', y='comparison_2',
                            hue='comparison_type', palette=['#7dc754', '#faa719', '#A680B8'],
                            linewidth=0.5, alpha=0.4)
    ax.plot([0, axis_lim], [0, axis_lim], 'red', linestyle='dashed', alpha = 0.6)
    ax.set(title='Genes with A2G identified from the\npairwise comparisons (T = ' + threshold + ')',
            xlabel='Comparison 1', ylabel='Comparison 2')
    plt.legend(title='Comparison type',
                labels=['_', 'WT_mCherry vs. WT_IMP2', 'WT_mCherry vs. mCherry_IMP2', 'WT_IMP2 vs. mCherry_IMP2'])
    plt.title('Genes with A2G identified from the\npairwise comparisons (T = ' + threshold + ')', 
                fontweight='bold')
    plt.gcf().set_size_inches(4, 4)
    fig.savefig(analyzed_result_path.joinpath('scatterplot_' + threshold + 'COPILOT.png'), dpi=400, bbox_inches='tight')
    plt.close(fig)
            
            
if __name__ == '__main__':
    workdir = sys.argv[1]
    threshold = sys.argv[2]
    region = sys.argv[3]

    if region != "":
        collapse_regions = [region] 
    else: 
        collapse_regions = ['site', 'window', 'CDS_UTR', 'transcript']
        
    for region in collapse_regions:                                                                    
        HTRIBE_result_path = workdir + '/' + region + '/OR/OR/result'
        comparison_type = ['WT_mCherry', 'WT_IMP2', 'mCherry_IMP2']
        comparison_pair = list(combinations(comparison_type, 2))

        # Create a folder storing results
        analyzed_result_path = Path(HTRIBE_result_path).joinpath('analyzed_result')
        create_new_folder(analyzed_result_path)
        genes_by_comparison_type = gather_genes_by_comparison_type(analyzed_result_path=analyzed_result_path, 
                                                                   HTRIBE_result_path=HTRIBE_result_path, 
                                                                   threshold=threshold)

        # SECTION 1: OVERLAPS BETWEEN PAIR COMPARISONS
        # Generate Venn diagram to compare the overlap between the detected genes for different comparisons
        generate_venn_diagram(genes_by_comparison_type_=genes_by_comparison_type,
                              analyzed_result_path_=analyzed_result_path, 
                              threshold_=threshold)
        # Report significant genes
        report_sig_genes(sig_gene_filename=analyzed_result_path.joinpath('sig_genes' + region + '.txt'), 
                         threshold=threshold, 
                         genes_by_comparison_type=genes_by_comparison_type, 
                         comparison_pair=comparison_pair)
                
        # SECTION 2: DISTANCE BETWEEN PAIR COMPARISONS
        # Report Jaccard distance based on shared genes
        main_output_filename = analyzed_result_path.joinpath('output.txt')
        write_pairwise_jaccard_distance(main_output_filename=main_output_filename, 
                                        threshold=threshold, 
                                        genes_by_comparison_type=genes_by_comparison_type, 
                                        comparison_pair=comparison_pair)
            
        # SECTION 3: COMPARE CURRENT RESULTS TO JEETS RESULTS
        compare_with_jeet_result(Jeet_result_path=workdir + '/Jeet_result/', 
                                 threshold=threshold, 
                                 genes_by_comparison_type=genes_by_comparison_type, 
                                 comparison_type=comparison_type, 
                                 main_output_filename=main_output_filename)
                
        # SECTION 4: PLOT COMMON GENES BETWEEN PAIR COMPARISONS
        coords_by_no_editing_sites_filename = analyzed_result_path.joinpath('coords_by_no_editing_sites' + threshold + '.pickle')
        coords_by_no_editing_sites = generate_coords_by_no_editing_sites(coords_by_no_editing_sites_filename=coords_by_no_editing_sites_filename, 
                                                                        comparison_pair=comparison_pair, 
                                                                        genes_by_comparison_type=genes_by_comparison_type)
        plot_scatterplots(coords_by_no_editing_sites=coords_by_no_editing_sites, 
                          analyzed_result_path=analyzed_result_path, 
                          threshold=threshold)
