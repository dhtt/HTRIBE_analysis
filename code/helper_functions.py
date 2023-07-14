import os 
import pandas as pd
import pandas as pd
import random
import numpy as np


def gather_read_counts(result_dir: str, normalize = True):
    '''This function reads the quantified expression from transcripts and return a data frame of concatenated results

    Args:
        result_dir (str): absolute path to the quantified transcripts result folder
        normalize (bool): whether the expression should be normalized by Deseq2 scheme or not
    '''
    all_files =['/'.join([result_dir, 'ctrl_DS' + str(file_no), 'quant.sf']) for file_no in range(1, 10, 1)]
    
    df_list = []
    for file in all_files:
        df = pd.read_csv(file, sep='\t', header=0, index_col=0)
        df_list.append(df[['NumReads']])
    count_matrix = pd.concat(df_list, axis=1)
    count_matrix.columns = ['_'.join([i, str(j)]) for i in ['WT', 'mCherry', 'IMP2'] for j in range(1, 4, 1)]
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
        os.mkdir(path_)
        
        

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