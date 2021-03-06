import numpy as np
import os
import json

# Transforms sequence in binary except "P" where "1" represents contiguous strong hydrophobic amino acids
# input: string sequence of AA
# output: string sequence in binary (except "P")
# P is Proline and is considered a "breaker"; therefore it is not binarized
def binarization(seq):
    # list of contiguous strong hydrophobic amino acids
    list_aa = ['V', 'I', 'L', 'F', 'M', 'Y', 'W', 'v', 'i', 'l', 'f', 'm', 'y', 'w']
    for i in range(len(seq)):
        if seq[i] in list_aa:
            yield '1'
        elif seq[i] in ['P', 'p']:
            yield 'P'
        else:
            yield '0'

# Find all the positions of an element in a sequence
# input: string sequence in binary (except "P"), string element to find
# output: generator of integer
def find_x(seq, x):
    for ind,nt in enumerate(seq):
        if nt == x:
            yield ind

# Delimitate all the HCs and find all the "1"-positions in them
# input: string sequence in binary (except the "P")
# output: generator of list
def find_HCs(seq):
    pos_1 = list(find_x(seq, '1'))
    pos_P = list(find_x(seq, 'P'))
    imp_pos = sorted(pos_1+pos_P)
    hc = []
    for pos in imp_pos:
        if pos not in pos_P and hc == []:
            hc.append(pos)
            continue
        elif pos not in pos_P and pos-hc[-1]<5:
            hc.append(pos)
            continue
        if len(hc)>1:
            yield hc
        if pos not in pos_P:
            hc = [pos]
        else:
            hc = []
    if len(hc)>1:
        yield hc

# Binarize all HCs
# input: list of "1"-position lists
# output: generator of HC sequence in binary
def make_hc(positions):
    for hc in positions:
        yield ''.join(['1' if i in hc else '0' for i in range(hc[0], hc[-1]+1)])

# Keep only the extremities of HCs (boundaries)
# input: list of position lists
# output: generator of start and end positions
def keep_limits(positions):
    for pos in positions:
        yield pos[0]
        yield pos[-1]

# Put back "." and "-" in the sequence
# input: list of HC limits, list of ponctuation positions
# output: list of adjusted HC limits
def put_back_ponctuation(positions, list_ponctuations):
    for pos in positions:
        for ponct in list_ponctuations:
            if pos >= ponct:
                pos+=1
            else:
                break
        yield pos

# Create a dictionary associating HC positions with the HC
# input: list of HC limits and list of HCs
# output: dictionary of positions tuples (start, end) associated to their HCs
def make_dict(positions, HCs):
    final_dict = {}
    for i in range(int(len(positions)/2)):
        final_dict[(positions[i*2], positions[i*2+1])] = HCs[i]
    return final_dict

# Find the different HCs in a sequence and give their associated positions
# input: sequence of AA
# output: dictionary of positions tuples (start, end) associated to their HCs
def binary_coding(seq):
    # remove "." and "-" from the sequence but save their position
    list_ponctuations = [ind for ind,nucl in enumerate(seq) if nucl in ['.', '-']]
    seq = seq.replace('.', '').replace('-', '')
    seq = ''.join(list(binarization(seq)))
    all_positions = list(find_HCs(seq))
    HCs = list(make_hc(all_positions))
    positions = list(keep_limits(all_positions))
    positions = list(put_back_ponctuation(positions, list_ponctuations))
    return make_dict(positions, HCs)

# Get the Pfam accession numbers
# input: file containing the soluble domains of known 3D structures + Pfam accession on last column
# output: list of Pfam accession numbers
def get_PF(file):
    f = open(file, 'r')
    list_PF = set()
    for line in f.readlines():
        list_PF.add(line.split('\t')[-1].strip())
    f.close()
    return list_PF

# Read the alignments of soluble domains in a database
# input: file containing alignments and file containing soluble domains
# output: dictionary associating superfamily to their alignment
def read_data(file_alignments, file_PF):
    list_PF = get_PF(file_PF)
    f = open(file_alignments, 'r', encoding = 'ISO-8859-1')
    data = {}
    data_i = []
    count=1
    for line in f.readlines():
        if line[0]=='/':
            if superfamily in list_PF:
                data[superfamily]=data_i
                data_i=[]
            count=1
        elif line[0]!='#' and superfamily in list_PF:
            data_i.append(line.split()[1])
        elif count==3:
            superfamily = line.split()[2].split('.')[0]
            count+=1
        elif count<5:
            count+=1
    f.close()
    return data

# Count how many each HC patterns are present in the soluble domains and update the count
# input: dictionary of positions tuples (start, end) associated to their HCs and dictionary associating each HC to its number of occurences
# output: dictionary associating each HC to its number of occurences
def count_HC(HCs, dico):
    for hc in HCs.values():
        dico.setdefault(int(hc, 2), 0)
        dico[int(hc, 2)]+=1
    return dico


# Get the Q-code of a HC
# input: integer of a HC
# output: Q-code string
def get_Q_codes(hc):
    hc=bin(hc)
    Q_code =''
    for i in range(len(hc)-1):
        if hc[i:i+2]=='11':
            Q_code+='V'
        elif hc[i:i+3]=='101':
            Q_code+='M'
        elif hc[i:i+4]=='1001':
            Q_code+='U'
        elif hc[i:i+5]=='10001':
            Q_code+='D'
    return Q_code

# Keep only the most numerous HCs
# input: dictionary associating HCs to its number of appearances and the number of HCs kept
# output: adjusted dictionary associating HCs to its number of occurences
def filtred_HCs(HCs, n):
    size = len(HCs)
    lim = int(np.quantile(list(HCs.values()), 1-n/size))
    HCs_filtred = {}
    for hc,count in HCs.items():
        if count>=lim:
            HCs_filtred[hc]=count
    return HCs_filtred

# Find most present HC patterns inside soluble domain alignments and optionally their Q-codes
# input: file containing alignments, file containing soluble domains, file for output, number of HCs kept and boolean which indicates if we need the Q-code
# output: file containing dictionary associating HCs with their number of occurences
def get_analyse(file_alignements, file_soluble_domains, output_file, nb_HC=500, Q_code=False):
    print('Preprocessing...')
    data = read_data(file_alignements, file_soluble_domains)
    print('done')
    print('Searching HCs...')
    all_HC = {}
    for _,alignment in data.items():
        for seq in alignment:
            HCs = binary_coding(seq)
            all_HC = count_HC(HCs, all_HC)
    print('done')
    
    print('Filtring...')
    all_HC = filtred_HCs(all_HC, nb_HC)
    all_HC = dict(sorted(all_HC.items(), key=lambda x:x[1], reverse=True))
    print('done')
    
    if Q_code:
        HC_with_Q_code = {}
        for hc,effectif in all_HC.items():
            HC_with_Q_code[(hc, get_Q_codes(hc))] = effectif
        return HC_with_Q_code
    
    # save results
    output = open(output_file, 'w')
    json.dump(all_HC, output)
    output.close()