import argparse
import sys
import os
import subprocess
import multiprocessing
from pathlib import Path
import pandas as pd
import json
from collections import defaultdict

def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", 
        "--input_path", 
        type=str,
        help="Path to project files."
    )
    parser.add_argument(
        "-r1", 
        "--read1_extension", 
        type=str,
        help="File extension for read 1."
    )
    parser.add_argument(
        "-r2", 
        "--read2_extension", 
        type=str,
        help="File extension for read 2."
    )
    parser.add_argument(
        "-b", 
        "--barcode", 
        type=str,
        default="N"*16,
        help="Barcode for UMI tool"
    )
    parser.add_argument(
        "-gen", 
        "--genome_dir", 
        type=str,
        help="Folder where genome are indexed"
    )
    parser.add_argument(
        "--extra_umi_extract_options", 
        type=str,
        help="Extra arguments for UMI extraction"
    )
    parser.add_argument(
        "--extra_trimming_options", 
        type=str,
        help="Extra arguments for trimming"
    )
    parser.add_argument(
        "--extra_alignment_options", 
        type=str,
        help="Extra arguments for STAR alignment"
    )
    parser.add_argument(
        "--extra_umi_dedup_options", 
        type=str,
        help="Extra arguments for UMI deduplication"
    )
    return parser.parse_args(args)


def execute_generic_function(script_):
    # print(script_)
    subprocess.call(script_, shell=True)


def generate_command_lines(sample_id, input_path, output_path, read1_extension, read2_extension, analysis_step, **kwargs):
    output_path.mkdir(parents=True, exist_ok=True)
    
    forward_read = sample_id + read1_extension
    reverse_read = sample_id + read2_extension
    script = ""
    
    for parameter in kwargs:
        kwargs[parameter]
    
    try:
        all_files = os.listdir(input_path)
        assert (
            all([file in all_files for file in [forward_read, reverse_read]]) == True
            ), "Check if pair-end files exist for %s." %(sample_id)
        
    except Exception as e:
        print(e)
        
    else: 
        if analysis_step == "extract":    
            try: 
                assert ("barcode" in kwargs), "Check barcode!"
            except Exception as e:
                print(e)
            else:
                barcode = kwargs["barcode"]
                script = "umi_tools extract \
                        --bc-pattern=%s \
                        --stdin %s \
                        --stdout %s \
                        --read2-in %s \
                        --read2-out=%s" %(
                        barcode,
                        input_path.joinpath(forward_read), 
                        output_path.joinpath(forward_read),
                        input_path.joinpath(reverse_read),
                        output_path.joinpath(reverse_read))
                    
        if analysis_step == "trim":
            script = "cutadapt \
                --trim-n --match-read-wildcards -n 4 -u 15 \
                -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
                -a CCAAGTCT -a ATCTCGTATGCCGTCTTCTGCTTG \
                -a GGGGGGGGGG -a AAAAAAAAAA  -A TTTTTTTTTT \
	            -e 0.2 -j 0 --nextseq-trim 20 -m 40:40 \
                -o %s -p %s %s %s" %(
                output_path.joinpath(forward_read), 
                output_path.joinpath(reverse_read),
                input_path.joinpath(forward_read),
                input_path.joinpath(reverse_read))
                
        if analysis_step == "align":
            try: 
                assert ("genome_dir" in kwargs), "Check if genome is indexded!"
            except Exception as e:
                print(e)
            else:
                genome_dir = kwargs["genome_dir"]
                script = "$STAR \
                    --genomeDir %s \
                    --runThreadN 50 \
                    --readFilesIn %s %s \
                    --outFileNamePrefix %s" %(
                    genome_dir, 
                    input_path.joinpath(forward_read), input_path.joinpath(reverse_read),
                    output_path.joinpath(sample_id + "_")
                )
                
        if analysis_step == "dedup":
            script = "umi_tools dedup \
                --stdin=%s \
	            --in-sam --out-sam > %s"%(
                input_path.joinpath(sample_id + "_"), #TODO
                output_path.joinpath(sample_id + ".sam")
            )
             
        return script


def execute_workflow(args=None):
    args = parse_args(args)
    print(args)
    
    input_path = Path(args.input_path)
    raw_files_path = input_path.joinpath('raw_files')
    processed_files_path = input_path.joinpath('processed')
    read1_extension = args.read1_extension
    read2_extension = args.read2_extension
    barcode = args.barcode
    genome_dir = input_path.joinpath(args.genome_dir)


    try:
        raw_files = os.listdir(raw_files_path)
        assert (len(raw_files) != 0), "Check if raw fastq files exist."
        
        for result_folder in ['extracted', 'extracted.trimmed', 'extracted.trimmed.aligned', 'extracted.trimmed.aligned.dedup']:
            new_folder = processed_files_path.joinpath(result_folder)
            new_folder.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        print(e)
        
    
    print ("========== Initialized workflow ==========")
    sample_ids = [file.replace(read1_extension, '').replace(read2_extension, '') for file in os.listdir(raw_files_path) \
        if file.endswith('fastq')]
    sample_ids = sorted(list(set(sample_ids)))
    print(sample_ids)    
    
    
    input_path = raw_files_path
    output_path = processed_files_path.joinpath('extracted')
    output_path.mkdir(parents=True, exist_ok=True)
    command_list = defaultdict(list)
    
    for sample_id in sample_ids:
        command_line = generate_command_lines(sample_id=sample_id, input_path=input_path, output_path=output_path, 
             read1_extension=read1_extension, read2_extension=read2_extension, analysis_step="extract", 
             barcode=barcode)
        command_list['extract'].append(command_line)
        
    with multiprocessing.Pool() as pool:
        pool.map(execute_generic_function, command_list['extract'])
        
        # command_line = generate_command_lines(sample_id=sample_id, input_path=input_path, output_path=output_path, 
        #      read1_extension=read1_extension, read2_extension=read2_extension, analysis_step="trim")
        
        # command_line = generate_command_lines(sample_id=sample_id, input_path=input_path, output_path=output_path, 
        #      read1_extension=read1_extension, read2_extension=read2_extension, analysis_step="align", 
        #      genome_dir=genome_dir)
        
        # command_line = generate_command_lines(sample_id=sample_id, input_path=input_path, output_path=output_path, 
        #      read1_extension=read1_extension, read2_extension=read2_extension, analysis_step="dedup")
        
    print("========== Finished ==========")


if __name__ == "__main__":
    sys.exit(execute_workflow())
