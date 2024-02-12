"""
This script is used to annotate A2G sites based on a configuration file.
It reads the configuration file, extracts the necessary arguments, and performs various operations on the A2G sites.
The script utilizes helper functions from the 'helper_functions' module.
"""

from collections import defaultdict
import subprocess
import argparse
import configparser
import os
from pathlib import Path
from helper_functions import *

def main():
    """
    The main function of the script.
    It parses the command-line arguments, reads the configuration file, and performs the A2G site annotation.
    """
    # Create an ArgumentParser instance
    parser = argparse.ArgumentParser(description='Call bedtools intersect with a reference path')
    parser.add_argument('config_file', type=str)
    args = parser.parse_args()
    
    # Read configuration from the specified config file
    config = configparser.ConfigParser()
    config.read(args.config_file)
    
    # Extract arguments from the configuration file
    span_type = config.get('METHOD', 'span_type')
    HTRIBE_result_path = config.get('PATHS', 'HTRIBE_result_path')
    reference_genome = config.get('PATHS', 'reference_genome')
    
    collapsing_modes = ['AND', 'OR']
    thresholds = ['', '_1%', '_5%']
    exp_groups = {'wt': [1, 2, 3], 'mcherry': [4, 5, 6], 'imp2': [7, 8, 9]}
    
    # Define paths
    os.chdir(HTRIBE_result_path)
    annotated_A2G_path = '/'.join(['all_span', span_type])
    
    # Perform collapsing and annotation operations on A2G sites
    for collapsing_mode in collapsing_modes:    
        collapsed_replicates_1_path = Path(annotated_A2G_path).joinpath(collapsing_mode)
        collapsed_replicates_2_path = Path(collapsed_replicates_1_path).joinpath(collapsing_mode)
        for comparison_type in exp_groups:
            summarized_result_path = Path(collapsed_replicates_2_path).joinpath('comparisons').joinpath(comparison_type)
            create_new_folder(summarized_result_path)

    # Collapse replicates
    all_files = os.listdir(HTRIBE_result_path)
    for file in all_files:
        if file.endswith('bedgraph'):
            remove_INTRON(file)
            if span_type in ['CDS_UTR', 'window', 'transcript']:
                annotated_bedgraph = '/'.join([annotated_A2G_path, file])
                if span_type == 'window':
                    bedtools_slop(file, reference_genome, 10, output_dir=annotated_A2G_path)
                else:
                    bedtools_intersect(reference_genome, file, annotated_bedgraph, '-wa')
                    bedtools_sort(annotated_bedgraph)
                    bedtools_groupby(annotated_bedgraph, '1-6', '7', 'distinct')
            elif span_type in ['site']:
                bedtools_sort(file, output_dir=annotated_A2G_path)
    
    all_annotated_files = os.listdir(annotated_A2G_path)
    for threshold in thresholds:
        print('*'*20 + threshold + '*'*20)
        
        intra_rep_filenames = defaultdict()
        for i, exp_group_1 in enumerate(exp_groups.keys()):
            for j, exp_group_2 in enumerate(exp_groups.keys()):
                # Compare only wt to mcherry, wt to imp2 and mcherry to imp2, not other way around
                if i < j: 
                    comparison_pair = exp_group_1 + '_' + exp_group_2
                    intra_rep_filenames[comparison_pair] = dict()
                    exp_group_1_idx = exp_groups[exp_group_1]
                    exp_group_2_idx = exp_groups[exp_group_2]
                    for idx_1 in exp_group_1_idx:
                        rep_ids = ['_' + str(idx_1) + '_' + str(idx_2) + '_' for idx_2 in exp_group_2_idx]
                        rep_filenames = [file for file in all_annotated_files if any(rep_id in file for rep_id in rep_ids) & file.endswith(threshold + '.bedgraph')]
                        id = str(idx_1) + '_' + ''.join([str(idx_2) for idx_2 in exp_group_2_idx])
                        intra_rep_filenames[comparison_pair][id] = sorted(rep_filenames)

        if span_type in ['CDS_UTR', 'transcript']:
            for comparison_pair, replicates in intra_rep_filenames.items():
                comparison_filenames = [item for items in replicates.values() for item in items]
                output_path = comparison_pair + '_A2G' + threshold + '.bedgraph.cat'
                concat_bedfiles(comparison_filenames, output_path)
        
        for collapsing_mode in collapsing_modes:   
            for comparison_pair, replicates in intra_rep_filenames.items():
                intra_group_collapse_output_dir = '/'.join(['all_span', span_type])
                for replicate, rep_filenames in replicates.items():
                    input_names = ['/'.join([annotated_A2G_path, rep_filename]) for rep_filename in rep_filenames]
                    output_name = comparison_pair + '_' + replicate + '_A2G' + threshold + '.bedgraph'
                    if len(rep_filenames) > 0:
                        find_dup_regions(input_names, intra_group_collapse_output_dir, output_name, collapsing_mode, span_type)            
                    
                inter_group_collapse_input_dir = intra_group_collapse_output_dir
                inter_rep_filenames = ['/'.join([inter_group_collapse_input_dir, file]) 
                                        for file in os.listdir(inter_group_collapse_input_dir) 
                                        if (comparison_pair in file) & file.endswith(threshold + '.bedgraph')]
                inter_group_collapse_output_dir = '/'.join([intra_group_collapse_output_dir, collapsing_mode])
                output_name = comparison_pair + '_A2G' + threshold + '.bedgraph'
                if len(inter_rep_filenames) > 0:
                    find_dup_regions(inter_rep_filenames, inter_group_collapse_output_dir, output_name, collapsing_mode, span_type)
        
            final_result_path = '/'.join([inter_group_collapse_output_dir, collapsing_mode])
            for bedgraph_file in os.listdir(final_result_path):
                if bedgraph_file.endswith('.bedgraph'):
                    if span_type not in ['site', 'window']:
                        result_input_path = bedgraph_file + '.site'
                        intersect_input_path = '/'.join([final_result_path, bedgraph_file])
                        intersect_output_path = '/'.join([final_result_path, result_input_path])
                        bedtools_intersect(bedgraph_file + '.cat', intersect_input_path, intersect_output_path, '-wa')
                    else:
                        result_input_path = bedgraph_file
                    get_HYPERTRIBE_result(result_input_path, final_result_path)


if __name__=='__main__':
    main()