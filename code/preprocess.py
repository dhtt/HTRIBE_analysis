import argparse
import sys
import os
import subprocess
from pathlib import Path
import pandas as pd
import json


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
        assert (len(raw_files) != 0), "Check if raw fastq files exist"
        
        for result_folder in ['extracted', 'extracted.trimmed', 'extracted.trimmed.aligned', 'extracted.trimmed.aligned.dedup']:
            new_folder = processed_files_path.joinpath(result_folder)
            new_folder.mkdir(parents=True, exist_ok=True)
            
    except Exception as e:
        print(e)
        
    
    
    
    # print ("========== Initialized workflow ==========")
    
    # sample_ids = [file.replace(read1_extension, '').replace(read2_extension, '') for file in os.listdir(raw_files_path)]
    # sample_ids = sorted(list(set(sample_ids)))
    # print(sample_ids)
    
    
    # print("Extract UMIs")
    # for sample_id in sample_ids:
    #     input_path = raw_files_path
    #     output_path = processed_files_path.joinpath('extracted')
    #     output_path.mkdir(parents=True, exist_ok=True)
        
    #     forward_read = sample_id + read1_extension
    #     reverse_read = sample_id + read2_extension
    #     subprocess.call(
    #         "umi_tools extract \
    #             --bc-pattern = %s \
    #             --stdin        %s \
    #             --stdout       %s \
    #             --read2-in     %s \
    #             --read2-out =  %s"
    #         %(
    #             barcode,
    #             input_path.joinpath(forward_read), 
    #             output_path.joinpath(forward_read),
    #             input_path.joinpath(reverse_read),
    #             output_path.joinpath(reverse_read)
    #         ),
    #         shell=True)
        
    # print("Trimming using cutadapt")
    # for sample_id in sample_ids:
    #     input_path = working_dir.joinpath('extracted')
    #     output_path = processed_files_path.joinpath('extracted.trimmed')
    #     output_path.mkdir(parents=True, exist_ok=True)
        
    #     forward_read = sample_id + read1_extension
    #     reverse_read = sample_id + read2_extension
    #     subprocess.call(
    #         "cutadapt \
    #             --trim-n --match-read-wildcards -n 4 -u 15 \
    #             -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    #             -a CCAAGTCT -a ATCTCGTATGCCGTCTTCTGCTTG \
    #             -a GGGGGGGGGG -a AAAAAAAAAA  -A TTTTTTTTTT \
	#             -e 0.2 -j 0 --nextseq-trim 20 -m 40:40 \
    #             -o %s \
    #             -p %s \
    #             %s %s"
    #         %(
    #             output_path.joinpath(forward_read), 
    #             output_path.joinpath(reverse_read),
    #             input_path.joinpath(forward_read),
    #             input_path.joinpath(reverse_read)
    #         ),
    #         shell=True)
        
    # print("STAR alignment")
    # for sample_id in sample_ids:
    #     input_path = working_dir.joinpath('extracted.trimmed')
    #     output_path = processed_files_path.joinpath('extracted.trimmed.aligned')
    #     output_path.mkdir(parents=True, exist_ok=True)
        
    #     forward_read = sample_id + read1_extension
    #     reverse_read = sample_id + read2_extension
    #     subprocess.call(
    #         "$STAR \
    #             --genomeDir %s \
    #             --runThreadN 50 \
    #             --readFilesIn %s %s \
    #             --outFileNamePrefix %s"
    #         %(
    #             genome_dir, 
    #             input_path.joinpath(forward_read), input_path.joinpath(reverse_read),
    #             output_path.joinpath(sample_id + "_")
    #         ),
    #         shell=True)
        
        
    # print("UMI Deduplication")
    # for sample_id in sample_ids:
    #     input_path = working_dir.joinpath('extracted.trimmed.aligned')
    #     output_path = processed_files_path.joinpath('extracted.trimmed.aligned.dedup')
    #     output_path.mkdir(parents=True, exist_ok=True)
        
    #     forward_read = sample_id + read1_extension
    #     reverse_read = sample_id + read2_extension
    #     subprocess.call(
    #         "umi_tools dedup \
    #             --stdin = %s \
	#             --in-sam --out-sam > %s ctrl_DS1.extract.trimmed_Aligned.out.dedup.sam"
    #         %(
    #             input_path.joinpath(sample_id + "_"), 
    #             output_path.joinpath(sample_id + ".sam")
    #         ),
    #         shell=True)
        
    print("========== Finished ==========")


if __name__ == "__main__":
    sys.exit(execute_workflow())
