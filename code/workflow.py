import os
import collections
import configparser
import argparse

import preprocess 

def get_config_section(config_ : configparser.RawConfigParser()):
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
        help="Path to experiment config file."
    )
    return parser.parse_args(args)


def run_DEU_workflow():
    preprocess.execute_workflow(args = [
        '--raw_files_path', raw_files_path,
        '--output_path', output_path,
        '--read1_extension', read1_extension, 
        '--read2_extension', read2_extension, 
        '--genome', genome,
        '--refgen_path', refgen_path,
        '--control_id', control_id,
        '--treatment_ids', treatment_ids,
        '--mrnaseq_options', mrnaseq_options
    ])


def run_DHM_workflow():
    DHM_workflow.execute_workflow(args = [
        '--chip_seq_path', chipseq_path,
        '--read1_extension', read1_extension_dhm, 
        '--read2_extension', read2_extension_dhm, 
        '--genome', genome,
        '--control_id', control_id,
        '--treatment_ids', treatment_ids,
        '--antibody', antibody,
        '--replicate_index', replicate_index,
        '--group_name_index', group_name_index,
        '--chipseq_options', chipseq_options,
        '--refgen_flank_path', refgen_flank_path
    ])

def get_DEU_DHM_correlation():
    correlate_DEU_DHM.execute_workflow(args = [
        '--mrna_seq_path', rnaseq_path,
        '--chip_seq_path', chipseq_path,
        '--refgen_path', refgen_path,
        '--refgen_flank_path', refgen_flank_path,
        '--control_id', control_id,
        '--treatment_ids', treatment_ids
    ])


if __name__ == "__main__":
    # Set working directory & parse arguments
    os.chdir('/home/dhthutrang/episplice-pipeline')
    args = parse_args(args=None)
    config_path = args.config_path
    run_option = args.run_option

    # Read config file
    config = configparser.RawConfigParser()
    config.read(config_path)
    config_dict = get_config_section(config)
    
    # Define variables from configurations
    input_path = config_dict['OVERALL_CONFIG']['input_path']
    genome = config_dict['OVERALL_CONFIG']['genome']

    strandedness = config_dict['mRNA_seq_CONFIG']['strandedness']
    read1_extension_deu = config_dict['mRNA_seq_CONFIG']['read1_extension']
    read2_extension_deu = config_dict['mRNA_seq_CONFIG']['read2_extension']
    refgen_path = config_dict['mRNA_seq_CONFIG']['refgen_path']
    control_id = config_dict['mRNA_seq_CONFIG']['control_id']
    treatment_ids = config_dict['mRNA_seq_CONFIG']['treatment_ids']
    n_cores = config_dict['mRNA_seq_CONFIG']['n_cores']
    mrnaseq_options = config_dict['mRNA_seq_CONFIG']['mrnaseq_options']


    read1_extension_dhm = config_dict['ChIP_seq_CONFIG']['read1_extension']
    read2_extension_dhm = config_dict['ChIP_seq_CONFIG']['read2_extension']
    antibody = config_dict['ChIP_seq_CONFIG']['antibody']
    replicate_index = config_dict['ChIP_seq_CONFIG']['replicate_index']
    group_name_index = config_dict['ChIP_seq_CONFIG']['group_name_index']
    chipseq_options = config_dict['ChIP_seq_CONFIG']['chipseq_options']
    refgen_flank_path = config_dict['ChIP_seq_CONFIG']['refgen_flank_path']

    rnaseq_path = '/'.join([input_path, 'mRNA_seq'])
    chipseq_path = '/'.join([input_path, 'ChIP_seq'])


    if run_option == "both":
        run_DEU_workflow()
        run_DHM_workflow()
    elif run_option == "deu":
        run_DEU_workflow()
    elif run_option == "dhm":
        run_DHM_workflow()
    elif run_option == "correlation":
        get_DEU_DHM_correlation()