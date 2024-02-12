import os
import collections
import configparser
import argparse

import preprocess 

def get_config_section(config_):
    if not hasattr(get_config_section, 'section_dict'):
        get_config_section.section_dict = collections.defaultdict()
        
        for section in config_.sections():
            get_config_section.section_dict[section] = dict(config_.items(section))
    
    return get_config_section.section_dict


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-cfg", 
        "--config_path", 
        type=str,
        default="/home/dhthutrang/TRIBE/mRNA_seq/HTRIBE.cfg",
        help="Path to experiment config file."
    )
    return parser.parse_args(args)


def run_preprocessing():
    preprocess.execute_workflow(args = [
        '--input_path', input_path,
        '--read1_extension', read1_extension, 
        '--read2_extension', read2_extension, 
        '--genome_dir', genome_dir,
        '--barcode', barcode,
        '--extra_umi_extract_options', extra_umi_extract_options,
        '--extra_trimming_options', extra_trimming_options,
        '--extra_alignment_options', extra_alignment_options,
        '--extra_umi_dedup_options', extra_umi_dedup_options
    ])
    

if __name__ == "__main__":
    # Set working directory & parse arguments
    args = parse_args(args=None)
    config_path = args.config_path

    # Read config file
    config = configparser.RawConfigParser()
    config.read(config_path)
    config_dict = get_config_section(config)
    
    # Define variables from configurations
    input_path = config_dict['general_CONFIG']['input_path']
    os.chdir(input_path)

    read1_extension           = config_dict['preprocessing_CONFIG']['read1_extension']
    read2_extension           = config_dict['preprocessing_CONFIG']['read2_extension']
    genome_dir                = config_dict['preprocessing_CONFIG']['genome_dir']
    barcode                   = config_dict['preprocessing_CONFIG']['barcode']
    extra_umi_extract_options = config_dict['preprocessing_CONFIG']['extra_umi_extract_options']
    extra_trimming_options    = config_dict['preprocessing_CONFIG']['extra_trimming_options']
    extra_alignment_options   = config_dict['preprocessing_CONFIG']['extra_alignment_options']
    extra_umi_dedup_options   = config_dict['preprocessing_CONFIG']['extra_umi_dedup_options']

    run_preprocessing()