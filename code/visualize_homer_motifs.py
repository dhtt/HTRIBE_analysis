import os
import argparse
import pandas as pd
from collections import defaultdict
import re 
import matplotlib.pyplot as plt
from sklearn.cluster import AffinityPropagation
from scipy.spatial.distance import hamming, pdist
from itertools import combinations, starmap
from functools import reduce
import numpy as np

def parse_motif_file(motif_file_path: str):
    with open(motif_file_path, 'r') as f:
        lines = f.readlines()
        f.close()
        
    folder_name = motif_file_path.split('/')[-2]
    print(folder_name)
    region = folder_name.split('.')[-1] if '.' in folder_name else 'all'
    comparison = '_'.join(folder_name.split('_')[:2]) 

    motif_info = dict()
    id, alt_name, pvalue, pos = [],[],[],[]
    
    for line in lines:
        if line.startswith('>'):
            motif = line.split('\t')
            id.append(motif[1])
            alt_name.append(motif[0].strip('>'))
            stats = motif[5]
            pvalue.append(stats.split(',')[-1].split(':')[-1])
            pos.append(motif[6])
    motif_info = {'id': id, 'alt_name': alt_name, 'pvalue': pvalue, 'pos': '-', 'region': region, 
                  'comparison': comparison} 
    return(motif_info)
         
def match_rule(char_):   
    rule_R = ['R', 'A', 'G']
    rule_D = ['D', 'A', 'G', 'U']
    rule_H = ['H', 'A', 'U', 'C']
    rule_Y = ['Y', 'C', 'U']
    rule_N = ['A', 'U', 'C', 'G']
    
    if char_ in rule_R:
        return('R')
    elif char_ in rule_H:
        return('H')
    elif char_ in rule_D:
        return('D')
    elif char_ in rule_Y:
        return('Y')
    elif char_ == '1':
        return('ACH')
    elif char_ == '2':
        return('AC')
    elif char_ == '3':
        return('AC')
    elif char_ == '4':
        return('AY')
    elif char_ == '5':
        return('UUU')
    elif char_ == 'u':
        return('U')
    else:
        return(char_)
    
def classify_motif(motif: str):
    motif = motif.replace('T', 'U')
    
    alias = '' 
    if ('ACH' in motif) | ('ACG' in motif): #['RRACH', 'RACH', 'DRACG', 'DRACH', 'URACH']
        motif_ = motif.replace('ACH', '1').replace('ACG', '2')
        for char in motif_:
            alias += match_rule(char)
    else:
        if ('AC' in motif) | ('AU' in motif) | ('AY' in motif): # ['RRAC', 'DRAY']
            motif_ = motif.replace('AC', '3').replace('AU', 'AY').replace('AY', '4')
            for char in motif_:
                alias += match_rule(char)
        else:
            if 'UUU' in motif:
                motif_ = motif.replace('UUU', '5')
                for char in motif_:
                    alias += match_rule(char)
            else:
                if re.search(r'U.U', motif):
                    motif_ = motif.replace('U', 'u')
                    alias = ['N' if char != 'u' else 'U' for char in motif_ ]
                else:
                    alias = ['N' for char in motif]
        
    return(''.join(alias))  

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Gather all homer motifs')   
    argparser.add_argument('--homer_motifs_path', type=str, help='Path to all homer motifs')
    args = argparser.parse_args()
    homer_motifs_path = args.homer_motifs_path
    
    all_motifs = []
    for root, dirs, files in os.walk(homer_motifs_path):
        if '1%' in root:            
            if 'homerMotifs.all.motifs' in files:
                motif_info = parse_motif_file("/".join([root, 'homerMotifs.all.motifs']))
                all_motifs.append(motif_info)
    all_motifs_df = pd.concat([pd.DataFrame.from_dict(df) for df in all_motifs]).reset_index(drop=True)
    all_motifs_df['motif_class'] = all_motifs_df['alt_name'].apply(lambda x: classify_motif(x))
    
    consensus = '|'.join(['RRACH', 'RRAC', 'RACH', 'DRACG', 'DRACH', 'DRAY', 'URACH', 'UNUNU', 'U.U'])
    all_motifs_df['match_consensus'] = all_motifs_df['motif_class'].apply(lambda x: True if re.search(consensus, x) else False)
    
    all_motifs_df.to_csv('/'.join([homer_motifs_path, 'all_homer_motifs.csv']), index=False, sep='\t')  
    
    
    
    
    
    # plt.close("all")
    # plt.figure(1)
    # plt.clf()

    # plt.title("Estimated number of clusters")
    # plt.savefig('/'.join([homer_motifs_path, 'motif_clustering.png']))
    