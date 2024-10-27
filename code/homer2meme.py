import argparse

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description='Convert homer motif file to meme motif file')   
    argparser.add_argument('--homer_motif_file', type=str, help='Homer motif file')
    args = argparser.parse_args()
    homer_motif_file = args.homer_motif_file
    meme_motif_file = '/'.join(['/'.join(homer_motif_file.split('/')[:-1]), homer_motif_file.split('/')[-1] + '.meme.motif'])
    print(meme_motif_file)
    
    with open(homer_motif_file, 'r') as f:
        lines = f.readlines()
        f.close()
        new_lines = []
        motif_length = 5
        alphabet = ''
        for line in lines:
            if line.startswith('>'):
                motif = line.split('\t')
                motif_id = motif[1]
                motif_alt_name = motif[0].strip('>')
                new_lines.append(f'\nMOTIF {motif_id} {motif_alt_name}\n')
                alphabet += motif_alt_name
                
                alphabet_length = 4
                source_sites = 20
                source_E_value = motif[5].split(',')[-1].split(':')[-1]
                new_lines.append(f'letter-probability matrix: alength= {alphabet_length} w= {motif_length} nsites= {source_sites} E= {source_E_value}\n')
                
                motif_length = 0
            else:
                new_lines.append(line)
                motif_length += 1
                
    with open(meme_motif_file, 'w') as f:
        rules = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'K': 'GT', 'M': 'AC', 'S': 'CG', 'W': 'AT', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
        alphabet_set = list(set([char for char in alphabet]))
        alphabet_set = '\n'.join([f'{char} = {rules[char]}' for char in alphabet_set if char in rules.keys()])
        
        f.write('MEME version 4\n\n')
        f.write(f'ALPHABET= ACGT\n\n')
        # f.write(f'ALPHABET "new_alphabet" DNA-like\n{alphabet_set}\nEND ALPHABET\n\n')
        f.write('strands: + -\n\n')
        f.write('Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n')
        f.writelines(new_lines)
        f.close()
        