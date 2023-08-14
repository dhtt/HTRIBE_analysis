import os
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path


def mkdir_p(path_to_create):
    """Create a directory

    Args:
        path_to_create (str): directory to be created
    """
    Path(path_to_create).mkdir(parents=True, exist_ok=True)


def rename_file(file, new_file):
    """Rename a file

    Args:
        file (str): old file name
        new_file (str): new file name
    """
    subprocess.call(
        "mv %s %s" % (file, new_file),
        shell=True
    )


def sort_and_clean_bed_file(bed_file):
    """Sort and clean bed file

    Args:
        bed_file (str): name of the bed file to be sorted
    """
    subprocess.call(
        "sort -k1,1 -k2,2n %s > %s"
        % (bed_file, bed_file + ".sorted"),
        shell=True
    )
    rename_file(bed_file + '.sorted', bed_file)


def bedtools_groupby(bed_file, group_by, grouped_columns, group_mode):
    """This method calls bedtools groupby

    Args:
        bed_file (str): name of bedgraph file for which the regions should be collapsed 
        group_by (str): see bedtools groupby
        grouped_columns (str): see bedtools groupby
        group_mode (str): see bedtools groupby 
    """
    command = "bedtools groupby -i %s -g %s -c %s -o %s > %s" % (
        bed_file, group_by, grouped_columns, group_mode, bed_file + '.grouped')
    print(command)
    subprocess.call(command, shell=True)
    rename_file(bed_file + '.grouped', bed_file)


def bedtools_multiinter(bed_files: [str], output_path: str, multiinter_extra_args: str):
    """This method calls bedtools multiinter

    Args:
        bed_files ([str]): list of bedgraph files for which the mutual regions should be determined 
        output_path (str): where output should be saved to
        multiinter_extra_args (str): extra arguments for bedtools multiinter. Defaults to ''.
    """
    bed_files = ' '.join(bed_files)
    command = "bedtools multiinter -i %s %s > %s" % (
        bed_files, multiinter_extra_args, output_path)
    print(command)
    subprocess.call(command, shell=True)


def bedtools_intersect(bed_file_1, bed_file_2, output_path,  intersect_extra_args: str):
    """This method calls bedtools intersect

    Args:
        bed_file_1 (str): file name 1
        bed_file_2 (str): file name 2
        output_path (str): where output should be saved to
        intersect_extra_args (str): extra arguments for bedtools intersect. Defaults to ''.
    """
    if isinstance(bed_file_2, list):
        bed_file_2 = ' '.join(bed_file_2)

    command = "bedtools intersect -a %s -b %s %s > %s" % (
        bed_file_1, bed_file_2, intersect_extra_args, output_path)
    print(command)
    subprocess.call(command, shell=True)


def find_dup_regions(bedgraphs, output_dir, output_name, collapse_mode, multiinter_extra_args='', intersect_extra_args=''):
    """Given a collapsing mode (AND or OR), this method will find the overlapping sites/regions between more than 2 (OR)
    or all (AND) replicates within the first experiment group in the comparison. The A2G sites/regions are filtered using 
    these mutual regions and are consensus sites/regions across replicates. This overlap process repeats for the second 
    experimental group. For example: 
    - Group 1 with ids [1, 2] and group 2 with ids [3, 4 ] -> All comparisons: 1_3, 1_4, 2_3, 2_4
    - First overlap:  1_34, 2_34
    - Second overlap: 12_34

    Args:
        bedgraphs (list[str]): a list of bedgraph files to be overlapped
        output_dir (str): where output should be saved to
        output_name (str): name of the resulting bedgraph
        collapse_mode (str): 'AND' or 'OR'
        multiinter_extra_args (str, optional): extra arguments for bedtools multiinter. Defaults to ''.
        intersect_extra_args (str, optional): extra arguments for bedtools intersect. Defaults to ''.
    """
    try:
        # Get the mutual regions between HTRIBE replicates that contain A2G sites
        multiinter_output_path = '/'.join([output_dir,
                                          output_name + '.multiinter'])
        bedtools_multiinter(
            bedgraphs, multiinter_output_path, multiinter_extra_args)

        # Intersect HTRIBE replicates for all mutual regions
        intersect_output_path = '/'.join([output_dir, output_name])
        bedtools_intersect(bedgraphs[0], bedgraphs[1:len(
            bedgraphs)], intersect_output_path, intersect_extra_args)
        bedtools_groupby(intersect_output_path, '1-24', '1', 'first')

        # Filter regions that are mutual in ALL (AND) or 2/3 (OR) replicates:
        annot_output_path = '/'.join([output_dir,
                                     output_name + '.' + collapse_mode])
        if collapse_mode == 'AND':
            filter_span = "grep '1,2,3' " + multiinter_output_path + " > " + annot_output_path
        elif collapse_mode == 'OR':
            filter_span = "grep '1,2\|1,3\|2,3' " + multiinter_output_path + " > " + annot_output_path
        print(filter_span)
        subprocess.run(filter_span, shell=True, check=True)

        # Group all replicates #1
        grouped_output_path = '/'.join([output_dir,
                                       collapse_mode, output_name])
        bedtools_intersect(intersect_output_path, annot_output_path,
                           grouped_output_path, intersect_extra_args)

    except subprocess.CalledProcessError as e:
        print(e.stderr)


def get_HYPERTRIBE_result(bedgraph, input_dir, only_exon: False):
    """This function get the summary of HYPERTRIBE result by calling the provided perl script. If only exon should
    be considered for the converted sites, filtering will be applied first to the bedgraph files containing A2G sites

    Args:
        bedgraph (str): name of bedgraph file with A2G sites to summarize. bedgraph contain no '.' and ends with '.bedgraph'
        input_dir (str): directory to the input bedgraph file 
        only_exon (False): if True return only sites found within exons
    """
    try:
        if only_exon:
            input_path = '/'.join([input_dir, bedgraph])
            output_dir = '/'.join([input_dir, 'EXON'])
            output_path = '/'.join([output_dir, bedgraph])
            mkdir_p(output_dir)

            filter_exon = "grep 'EXON' " + input_path + " > " + output_path
            print(filter_exon)
            subprocess.run(filter_exon, shell=True, check=True)

            input_path_exon_only = output_path
            output_path_exon_only = '/'.join([output_dir,
                                             bedgraph.split('.')[0] + '.xls'])
            command = "perl $HYPERTRIBE/summarize_results.pl %s > %s" % (
                input_path_exon_only, output_path_exon_only)
            subprocess.call(command, shell=True)
        else:
            input_path = '/'.join([input_dir, bedgraph])
            output_path = '/'.join([input_dir,
                                   bedgraph.split('.')[0] + '.xls'])
            command = "perl $HYPERTRIBE/summarize_results.pl %s > %s" % (
                input_path, output_path)
            subprocess.call(command, shell=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr)


def read_bedgraph(result_dir: str, threshold: str, file_extension: str = '.xls',
                  header: int = 0, index=None):
    '''This function reads HYPERTRIBE result files in xls format and return a data frame of concatenated results

    Args:
        result_dir (str): absolute path to the HYPERTRIBE result folder
        threshold (str): threshold of HYPERTRIBE analysis to be analyzed. This should be a string starting with _A2G.
            For example, '_A2G_1%'
    '''
    all_files = os.listdir(result_dir)
    df_list = []
    for file in all_files:
        if file.endswith(threshold + file_extension):
            comparison_type = '_'.join(file.split('_')[1:3])
            df = pd.read_csv(os.path.join(result_dir, file),
                             sep='\t', header=header, index_col=index)
            df['comparison_type'] = comparison_type
            df_list.append(df)
    df_list_concat = pd.concat(df_list)
    return(df_list_concat)


def gather_read_counts(result_dir: str, normalize=True):
    '''This function reads the quantified expression from transcripts and return a data frame of concatenated results

    Args:
        result_dir (str): absolute path to the quantified transcripts result folder
        normalize (bool): whether the expression should be normalized by Deseq2 scheme or not
    '''
    all_files = ['/'.join([result_dir, 'ctrl_DS' + str(file_no), 'quant.sf'])
                 for file_no in range(1, 10, 1)]

    df_list = []
    for file in all_files:
        df = pd.read_csv(file, sep='\t', header=0, index_col=0)
        df_list.append(df[['NumReads']])
    count_matrix = pd.concat(df_list, axis=1)
    count_matrix.columns = ['_'.join([i, str(j)]) for i in [
        'WT', 'mCherry', 'IMP2'] for j in range(1, 4, 1)]
    print(count_matrix)

    if normalize:
        print(count_matrix)
        normed_count_matrix = deseq_normalize(count_matrix)
        return(normed_count_matrix.transpose())
    else:
        return(count_matrix.transpose())


def create_new_folder(path_):
    """Create folder if not exists

    Args:
        path_ (Path): Path to directory
    """
    if os.path.exists(path_) == False:
        try:
            path_.mkdir(parents=True)
        except OSError as e:
            print(e.strerror)


def deseq_normalize(count_df: pd.DataFrame):
    """Normalize gene expression using deseq

    Args:
        count_df (pd.DataFrame): transcripts count matrix
    """
    count_df = count_df[count_df.sum(axis=1) != 0] + 0.01
    pseudo_reference = count_df.prod(axis=1)**(1/count_df.shape[1])
    print(pseudo_reference)
    ratios = count_df.div(pseudo_reference, axis=0)
    print(ratios)
    medians = ratios.median(axis=0)
    print(medians)
    normed_counts = count_df.div((medians), axis=1)
    print(normed_counts)
    return(normed_counts)


# data = pd.DataFrame(np.random.randint(2, size=(10,3)), columns=list('ABC'))
# print(data)
# deseq_normalize(data)
