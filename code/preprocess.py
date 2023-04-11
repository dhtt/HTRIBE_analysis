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
    print(script_)
    # subprocess.call(script_, shell=True)


def generate_command_lines(sample_id, input_path, output_path, read1_extension, read2_extension, analysis_step, **kwargs):
    
    forward_read = sample_id + read1_extension
    reverse_read = sample_id + read2_extension
    script = "undefined"
    
    for parameter in kwargs:
        kwargs[parameter]
    
    try:
        all_files = os.listdir(input_path)

        # Check if paired-end files exist
        if analysis_step != "dedup":
            assert (
                all([file in all_files for file in [forward_read, reverse_read]]) == True
                ), "Check if pair-end files exist for %s." %(sample_id)
        # else:
        #     align_file = sample_id + "_align.out.sam"
        #     print(align_file)
        #     print(all_files)
        #     assert (
        #         align_file in all_files == True
        #         ), "No alignment file for %s." %(sample_id)
    except Exception as e:
        print(e)
        
    else: 
        # Define file requirements for each analysis task and write scripts
        if analysis_step == "extract":    
            try:
                assert ("barcode" in kwargs), "UMI-extract requires custom barcode(s). Check them!"
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
                    --runThreadN 20 \
                    --readFilesIn %s %s \
                    --outFileNamePrefix %s" %(
                    genome_dir, 
                    input_path.joinpath(forward_read), input_path.joinpath(reverse_read),
                    output_path.joinpath(sample_id + "_")
                )
                
        if analysis_step == "dedup":
            script = "umi_tools dedup \
                --stdin=%s \
                --log=%s\
	            --in-sam --out-sam > %s"%(
                input_path.joinpath(sample_id + "_Align.out.sam"), 
                input_path.joinpath(sample_id + ".dedup.log"), 
                output_path.joinpath(sample_id + ".sam")
            )
             
        return script


def execute_workflow(args=None):
    args = parse_args(args)
    print(args)
    
    # Define variables by arguments from config file
    input_path = Path(args.input_path)
    raw_files_path = input_path.joinpath('raw_files')
    processed_files_path = input_path.joinpath('processed')
    read1_extension = args.read1_extension
    read2_extension = args.read2_extension
    barcode = args.barcode
    genome_dir = input_path.joinpath(args.genome_dir)

    # Check if raw files are in raw_files_path
    try:
        raw_files = os.listdir(raw_files_path)
        assert (len(raw_files) != 0), "Check if raw fastq files exist."
    except Exception as e:
        print(e)
        
    
    print ("========== Initialized workflow ==========")
    # Get sample IDS by removing extension indicating forward/reverse reads
    sample_ids = [file.replace(read1_extension, '').replace(read2_extension, '') for file in os.listdir(raw_files_path) \
        if file.endswith('.fastq') | file.endswith('.fastq.gz')]
    sample_ids = sorted(list(set(sample_ids)))
    print(sample_ids)    
    
    
    command_list = defaultdict(list)
    tasks = ['extract', 'trim', 'align', 'dedup']
    
    # For each task, create and execute bash script 
    for task_id, task in enumerate(tasks): 
        print(str(task_id) + " " + task)
        
        # Extract step takes data from raw_files_path and creates output to be stored in processed_files_path/extract
        if task == 'extract':
            input_path = raw_files_path
            output_path = processed_files_path.joinpath('extract')   
        else:
            # For later steps, output_path is named by completed tasks appended together and stored in 
            # processed_files_path. input_path is the path to output of the previous task 
            input_path = processed_files_path.joinpath('.'.join(tasks[:(task_id)]))
            output_path = processed_files_path.joinpath('.'.join(tasks[:task_id+1]))
    
        # If task has not been initiated
        if os.path.isdir(output_path) != False:
            # Create path to store outputs
            output_path.mkdir(parents=True, exist_ok=True)
        
            # Generate command lines for each sample
            for sample_id in sample_ids:
                command_line = generate_command_lines(sample_id=sample_id, input_path=input_path, output_path=output_path, 
                    read1_extension=read1_extension, read2_extension=read2_extension, analysis_step=task,
                    barcode=barcode, genome_dir=genome_dir)
                command_list[task].append(command_line)
                
            # Execute command lines for samples simultaneously 
            with multiprocessing.Pool() as pool:
                pool.map(execute_generic_function, command_list[task])
    print("========== Finished ==========")


if __name__ == "__main__":
    sys.exit(execute_workflow())
