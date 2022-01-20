import pandas as pd
import numpy as np
import json
import os
import matplotlib.pyplot as plt
from itertools import combinations
from IPython.display import clear_output
from math import *
from src.find_HC import *

# Binarizate a sequene specifically to have 1 at HC emplacements
# input: string sequence and dictionary containing positions of hcs on each sequence
# output: binary sequence
def light_HCs(seq, hcs):
    seq = [0]*len(seq)
    for pos in hcs:
        for aa in range(pos[0], pos[1]+1):
            seq[aa]=1
    return seq

def get_hc_from_dico(alignment):
    for HCs in alignment.values():
        for hcs in HCs:
            yield hcs

def find_conserved_regions(alignment, threshold):
    alignment_hc=list(get_hc_from_dico(alignment))
    coverages = np.mean(list(alignment_hc), axis=0)
    n = len(coverages)
    conserved_pos = []
    pos=[]
    for i in range(n):
        if coverages[i]>=threshold:
            pos.append(i)
        elif pos!=[]:
            conserved_pos.append(pos)
            pos=[]
    if pos!=[]:
        conserved_pos.append(pos)
    return conserved_pos

def find_region(pos, HCs):
    for pos_hc in HCs:
        if (pos[0]>=pos_hc[0] and pos[0]<=pos_hc[1]) or (pos[-1]>=pos_hc[0] and pos[-1]<=pos_hc[1]):
            yield HCs[pos_hc]

def find_corresponding_HCs(conserved_pos, seq, HCs):
    for pos in conserved_pos:
        yield list(find_region(pos, HCs[seq][0]))

def count_hc(substitutes):
    dict_substitutes = dict([(sub, 0) for sub in set(substitutes)])
    for sub in substitutes:
        dict_substitutes[sub]+=1
    return dict_substitutes

def add_count(matrix, substitutes):
    substitutes = count_hc(substitutes)
    for sub,count in substitutes.items():
        matrix[sub, sub]+=sum(range(count-1))
    for combi in combinations(substitutes.keys(), 2):
        matrix[combi[0], combi[1]]+=substitutes[combi[0]]+substitutes[combi[1]]
    return matrix

def update_matrix(matrix, conserved_pos, HCs, common_HCs):
    order = list(common_HCs.keys())
    n=len(conserved_pos)
    conserved_hcs = [[] for i in range(n)]
    w=len(HCs)
    for ind,seq in enumerate(HCs.keys()):
        for rep in range(len(HCs[seq])):
            corresponding_regions = list(find_corresponding_HCs(conserved_pos, seq, HCs))
            for i in range(n):
                for hc in corresponding_regions[i]:
                    if hc==None:
                        continue
                    if int(str(hc), 2) in common_HCs:
                        conserved_hcs[i].append(order.index(int(hc, 2)))
    for i in range(n):
        matrix = add_count(matrix, conserved_hcs[i])
    return matrix

def symetrize(mat):
    for i in range(len(mat)-1):
        for j in range(i+1, len(mat)):
            mat[j, i]=mat[i, j]
    return mat

def get_count_matrix(alignments_file, PF_file, HC_counts_file, threshold, output_file, pseudo_count=False):
    print('Preprocessing...')
    data = read_data(alignments_file, PF_file)
    file = open(HC_counts_file, 'r')
    common_HCs = json.load(file)
    file.close()
    print('done')
    print('Construction ...')
    nb_alignments = len(data)
    counter = 0
    n = len(common_HCs)
    HCs_conserved=0
    if pseudo_count:
        matrix = np.ones((n, n))
    else:
        matrix = np.zeros((n, n))
    for alignment in data.values():
        counter+=1
        print(counter, '/', nb_alignments)
        print(len(alignment), 'x', len(alignment[0]))
        HCs_alignment = {}
        alignment_HCs = {}
        for ind,seq in enumerate(alignment):
            HCs_seq = binary_coding(seq)
            if seq not in HCs_alignment:
                HCs_alignment[seq] = [light_HCs(seq, HCs_seq)]
                alignment_HCs[seq] = [HCs_seq]
            else: # some sequences may be in multiple copies in the alignment
                HCs_alignment[seq].append(light_HCs(seq, HCs_seq))
                alignment_HCs[seq].append(HCs_seq)
        conserved_pos = find_conserved_regions(HCs_alignment, threshold)
        HCs_conserved+=len(conserved_pos)
        matrix = update_matrix(matrix, conserved_pos, alignment_HCs, common_HCs)
    matrix = symetrize(matrix)
    print('done')
    # save results
    np.savetxt(output_file, matrix, delimiter=' ', fmt='%d')
    return len(data), HCs_conserved