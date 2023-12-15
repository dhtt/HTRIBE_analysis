#conda env python3.11
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


def remove_INTRON(input_path):
    command =  "grep 'EXON' " + input_path + " > " + input_path + '.exon'
    print(command)
    subprocess.call(command, shell=True)
    rename_file(input_path + '.exon', input_path)

    
def rename_file(file, new_file):
    """Rename a file

    Args:
        file (str): old file name
        new_file (str): new file name
    """
    command =  "mv %s %s" % (file, new_file)
    print(command)
    subprocess.call(command, shell=True)


def bedtools_slop(bed_file, reference_genome, window_length, output_dir=None, extra_sort_args=None):
    """Sort and clean bed file

    Args:
        bed_file (str): name of the bed file to be sorted
        output_dir (str): alternate directory to save bedfile
    """
    if output_dir is None:
        output_name = bed_file + ".slop"
    else:
        output_name = '/'.join([output_dir, bed_file + ".slop"])

    if extra_sort_args is None:
        extra_sort_args = ''
        
    command =  "bedtools slop -i %s -g %s -b %s > %s" % (bed_file, reference_genome, window_length, output_name)
    print(command)
    subprocess.call(command, shell=True)
    rename_file(output_name, output_name.removesuffix('.slop'))
    

def bedtools_sort(bed_file, output_dir=None, extra_sort_args=None, keep_name=True):
    """Sort and clean bed file

    Args:
        bed_file (str): name of the bed file to be sorted
        output_dir (str): alternate directory to save bedfile
    """
    if output_dir is None:
        output_name = bed_file + ".sorted"
    else:
        output_name = '/'.join([output_dir, bed_file + ".sorted"])

    if extra_sort_args is None:
        extra_sort_args = ''
        
    command =  "sort -k1,1 -k2,2n -k3,3n %s %s > %s" % (extra_sort_args, bed_file, output_name)
    print(command)
    subprocess.call(command, shell=True)
    if keep_name:
        rename_file(output_name, output_name.removesuffix('.sorted'))


# def bedtools_merge(bed_file, grouped_columns, group_mode): #
#     """This method calls bedtools groupby

#     Args:
#         bed_file (str): name of bedgraph file for which the regions should be collapsed 
#         group_by (str): see bedtools groupby
#         grouped_columns (str): see bedtools groupby
#         group_mode (str): see bedtools groupby 
#     """
#     command = "bedtools merge -i %s -c %s -o %s > %s" % (
#         bed_file, grouped_columns, group_mode, bed_file + '.grouped')
#     print(command)
#     subprocess.call(command, shell=True)
#     rename_file(bed_file + '.grouped', bed_file)


def bedtools_groupby(bed_file, group_by, retained_columns, group_mode, keep_name=True):
    """This method calls bedtools groupby

    Args:
        bed_file (str): name of bedgraph file for which the regions should be collapsed 
        group_by (str): see bedtools groupby
        grouped_columns (str): see bedtools groupby
        group_mode (str): see bedtools groupby 
    """
    
    
    if isinstance(group_by, str):
        grouped_cols = group_by
    elif isinstance(group_by, tuple):
        grouped_cols = str(group_by[0]) + '-' + str(group_by[1])
    
    if isinstance(retained_columns, str):
        retained_cols = retained_columns
    elif isinstance(retained_columns, tuple):
       retained_cols = ','.join([str(i) for i in range(retained_columns[0], retained_columns[1])])
    
    if len(group_mode.split(',')) == 1:
        grouping_mode = ','.join([group_mode] * len(retained_cols.split(',')))
    else:
        grouping_mode = group_mode
        
    command = "bedtools groupby -i %s -g %s -c %s -o %s > %s" % (
        bed_file, grouped_cols, retained_cols, grouping_mode, bed_file + '.grouped')
    print(command)
    subprocess.call(command, shell=True)
    if keep_name:
        rename_file(bed_file + '.grouped', bed_file)


def bedtools_merge(bed_file, grouped_columns, group_mode):
    """This method calls bedtools groupby

    Args:
        bed_file (str): name of bedgraph file for which the regions should be collapsed 
        group_by (str): see bedtools groupby
        grouped_columns (str): see bedtools groupby
        group_mode (str): see bedtools groupby 
    """
    
    if isinstance(grouped_columns, str):
        grouped_cols = grouped_columns
    elif isinstance(grouped_columns, tuple):
        grouped_cols = ','.join([str(i) for i in range(grouped_columns[0], grouped_columns[1])])
    
    if len(group_mode.split(',')) == 1:
        group_mode = ','.join([group_mode] * len(grouped_cols.split('')))
    else:
        group_mode = group_mode
    
    command = "bedtools merge -i %s -c %s -o %s > %s" % (
        bed_file, grouped_cols, group_mode, bed_file + '.grouped')
    print(command)
    subprocess.call(command, shell=True)
    rename_file(bed_file + '.grouped', bed_file)


def bedtools_multiinter(bed_files: [str], output_path: str, multiinter_extra_args: str = None):
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


def bedtools_intersect(bed_file_1, bed_file_2, output_path, intersect_extra_args: str = None):
    """This method calls bedtools intersect

    Args:
        bed_file_1 (str): file name 1
        bed_file_2 (str): file name 2
        output_path (str): where output should be saved to
        intersect_extra_args (str): extra arguments for bedtools intersect. Defaults to ''.
    """
    if isinstance(bed_file_2, list):
        bed_file_2 = ' '.join(bed_file_2)

    command = "bedtools intersect -a %s -b %s %s > %s" % (bed_file_1, bed_file_2, intersect_extra_args, output_path)
    print(command)
    subprocess.call(command, shell=True)


def concat_bedfiles(bed_files: [str], output_path: str, grouped_columns = (1, 9), retained_columns = (10, 25), group_mode = 'first',
                    bedtools_intersect_extra_args: str = None):
    """This method calls bedtools multiinter

    Args:
        bed_files ([str]): list of bedgraph files for which the mutual regions should be determined 
        output_path (str): where output should be saved to
        multiinter_extra_args (str): extra arguments for bedtools multiinter. Defaults to ''.
    """
    bed_files = ' '.join(bed_files)
    command = "cat %s > %s" % (bed_files, output_path)
    print(command)
    subprocess.call(command, shell=True)
    
    # Sort but using gene name too 
    bedtools_sort(output_path, extra_sort_args='-k6,6')
    
    grouped_cols = str(grouped_columns[0]) + '-' + str(grouped_columns[1])
    retained_cols = ','.join([str(i) for i in range(retained_columns[0], retained_columns[1])])
    grouping_mode = ','.join([group_mode] * (int(retained_columns[1]) - int(retained_columns[0])))
    bedtools_groupby(output_path, grouped_cols, retained_cols, grouping_mode)

def find_dup_regions(bedgraphs, output_dir, output_name, collapse_mode, span_type, multiinter_extra_args='', intersect_extra_args=''):
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
        # Intersect HTRIBE replicates for all mutual regions
        intersect_output_path = '/'.join([output_dir, output_name])
        bedtools_intersect(bedgraphs[0], bedgraphs[1:len(bedgraphs)], intersect_output_path, intersect_extra_args)
    
        if span_type in ['window', 'site']:
            bedtools_sort(intersect_output_path)
            retained_columns = ','.join([str(i) for i in range(10, 24)])
            bedtools_groupby(intersect_output_path, '1,4,5,6,7,8,9,24', '2,3,' + retained_columns, 'min,max,' + ','.join(['median']*14))
            
            sort_column = "awk 'BEGIN { FS = \"\t\"; OFS = \"\t\" } { print $1, $9, $10, %s, %s, $8}' %s > %s" %(
                ', '.join([ '$' + str(i) for i in range(2, 8)]), 
                ', '.join([ '$' + str(i) for i in range(11, 25)]),
                intersect_output_path, intersect_output_path + '.sorted_columns')
            print(sort_column)
            subprocess.call(sort_column, shell=True) 
            rename_file(intersect_output_path + '.sorted_columns', intersect_output_path)
        else:
            bedtools_sort(intersect_output_path)
            bedtools_merge(intersect_output_path, '4,5,6', 'distinct,first,first')

        # Get the mutual regions between HTRIBE replicates that contain A2G sites
        multiinter_output_path = '/'.join([output_dir, output_name + '.multiinter'])
        bedtools_multiinter(bedgraphs, multiinter_output_path, multiinter_extra_args)
        
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
        grouped_output_path = '/'.join([output_dir, collapse_mode, output_name])
        bedtools_intersect(intersect_output_path, annot_output_path, grouped_output_path, intersect_extra_args)


    except subprocess.CalledProcessError as e:
        print(e.stderr)


def get_HYPERTRIBE_result(bedgraph, input_dir):
    """This function get the summary of HYPERTRIBE result by calling the provided perl script. If only exon should
    be considered for the converted sites, filtering will be applied first to the bedgraph files containing A2G sites

    Args:
        bedgraph (str): name of bedgraph file with A2G sites to summarize. bedgraph contain no '.' and ends with '.bedgraph'
        input_dir (str): directory to the input bedgraph file 
        only_exon (False): if True return only sites found within exons
    """
    try:
        input_path = '/'.join([input_dir, bedgraph])
        output_dir = '/'.join([input_dir, 'result'])
        output_path = '/'.join([output_dir, bedgraph.split('.')[0] + '.xls'])
        mkdir_p(output_dir)

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
    print(all_files)
    df_list = []
    for file in all_files:
        if file.endswith(threshold + file_extension):
            comparison_type = '_'.join(file.split('_')[0:2])
            df = pd.read_csv(os.path.join(result_dir, file),
                             sep='\t', header=header, index_col=index)
            df['comparison_type'] = comparison_type
            df_list.append(df)
    df_list_concat = pd.concat(df_list)
    print(df_list_concat)
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
