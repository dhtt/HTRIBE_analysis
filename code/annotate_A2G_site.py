from collections import defaultdict
import subprocess
import argparse
import configparser
import os
from pathlib import Path
from helper_functions import create_new_folder, find_dup_regions, sort_and_clean_bed_file, get_HYPERTRIBE_result


def main():
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
    annotated_A2G_path = Path('all_span').joinpath(span_type)
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
            sort_and_clean_bed_file(file)
        
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
                        rep_ids = ['_' + str(idx_1) + '_' + str(idx_2) for idx_2 in exp_group_2_idx]
                        rep_filenames = [file for file in all_files if any(rep_id in file for rep_id in rep_ids) & file.endswith(threshold + '.bedgraph')]
                        id = str(idx_1) + '_' + ''.join([str(idx_2) for idx_2 in exp_group_2_idx])
                        intra_rep_filenames[comparison_pair][id] = sorted(rep_filenames)
        
        
        for collapsing_mode in collapsing_modes:   
            for comparison_pair, replicates in intra_rep_filenames.items():
                intra_group_collapse_output_dir = '/'.join(['all_span', span_type])
                for replicate, rep_filenames in replicates.items():
                    output_name = comparison_pair + '_' + replicate + '_A2G' + threshold + '.bedgraph'
                    # find_dup_regions(rep_filenames, intra_group_collapse_output_dir, output_name, collapsing_mode)
                    
                inter_group_collapse_input_dir = intra_group_collapse_output_dir
                inter_rep_filenames = ['/'.join([inter_group_collapse_input_dir, file]) 
                                        for file in os.listdir(inter_group_collapse_input_dir) 
                                        if (comparison_pair in file) & file.endswith(threshold + '.bedgraph')]
                inter_group_collapse_output_dir = '/'.join([intra_group_collapse_output_dir, collapsing_mode])
                output_name = comparison_pair + '_A2G' + threshold + '.bedgraph'
                # find_dup_regions(inter_rep_filenames, inter_group_collapse_output_dir, output_name, collapsing_mode)
        
            final_result_path = '/'.join([inter_group_collapse_output_dir, collapsing_mode])
            for bedgraph_file in os.listdir(final_result_path):
                if bedgraph_file.endswith('.bedgraph'):
                    get_HYPERTRIBE_result(bedgraph_file, final_result_path, only_exon=False)
                    get_HYPERTRIBE_result(bedgraph_file, final_result_path, only_exon=True)


if __name__=='__main__':
    main()