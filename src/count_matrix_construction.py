import numpy as np
import json
from itertools import combinations
from src.find_HC import *

###############################################################################################################################
# Binarize a sequene specifically to have 1 at HC emplacements
# input: string sequence and dictionary containing positions of hcs on each sequence
# output: binary sequence
def light_HCs(seq, hcs):
    seq = [0]*len(seq)
    for pos in hcs:
        for aa in range(pos[0], pos[1]+1):
            seq[aa]=1
    return seq

# Get all binarized sequences from a dictionary, especially those in repeted sequences
# input: dictionary associating sequence to a list containing the binarized sequence a number of times equal to the number of occurences
# output: generator of binarized sequences
def get_binarized_sequences(alignment):
    for HCs in alignment.values():
        for hcs in HCs:
            yield hcs

# Find conserved regions of an binarized sequences alignment (rate of HCs appearance on a position above a threshold)
# input: dictionary associating sequence to a list containing the binarized sequence a number of times equal to the number of occurences and a threshold between 0 and 1
# output: list of lists containing regions positions
def find_conserved_regions(alignment, threshold):
    alignment_hc=list(get_binarized_sequences(alignment))
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

###############################################################################################################################
# Find HCs corresponding to the positions of a conserved region 
# input: list of positions and dictionary containing HCs positions associated to the HC
# output: matching HC generator
def find_region(pos, HCs):
    for pos_hc in HCs:
        if (pos[0]>=pos_hc[0] and pos[0]<=pos_hc[1]) or (pos[-1]>=pos_hc[0] and pos[-1]<=pos_hc[1]):
            yield HCs[pos_hc]

# Find all HCs corresponding with conserved regions
# input: list of list containing all conserved positions, sequence of the alignment and dictionary associating a sequence with its HCs positions and HCs patterns
# output: generator of HCs matching with a conserved region
def find_corresponding_HCs(conserved_pos, seq, HCs):
    for pos in conserved_pos:
        yield list(find_region(pos, HCs[seq][0]))

# Get a the count of appearance of each HC in an list
# input: list of HCs
# output: dictionary of each HC associated with its number of appearances
def count_hc(substitutes):
    dict_substitutes = dict([(sub, 0) for sub in set(substitutes)])
    for sub in substitutes:
        dict_substitutes[sub]+=1
    return dict_substitutes

# Add the substitution counts to the matrix
# input: numpy matrix and list of HCs
# output: numpy matrix
def add_count(matrix, substitutes):
    substitutes = count_hc(substitutes)
    for sub,count in substitutes.items():
        matrix[sub, sub]+=sum(range(count-1))
    for combi in combinations(substitutes.keys(), 2):
        matrix[combi[0], combi[1]]+=substitutes[combi[0]]+substitutes[combi[1]]
    return matrix

# Update the matrix with the subsitution counts of an alignment
# input: numpy matrix, list of lists containing positions of a conserved position and list of the most common HCs
# output: numpy matrix
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

###############################################################################################################################
# Get a symmetric matrix
# input: triangular matrix
# output: symmetric matrix
def symmetrize(mat):
    for i in range(len(mat)-1):
        for j in range(i+1, len(mat)):
            mat[j, i]=mat[i, j]
    return mat

###############################################################################################################################
# Get the count matrix containing all the substitutions that appears in the different alignments on conserved regions
# input: path to file containing all the alignment, path to file containing all the soluble domains, path to file containing the result of an HC count (obtained by the function "get_analyse()" for example), the coverage threshold, path to the output file, boolean to have a pseudo-count or no
# output: file containing the matrix and number of conserved regions in all the alignments
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
    conserved_regions=0
    if pseudo_count:
        matrix = np.ones((n, n))
    else:
        matrix = np.zeros((n, n))
    for alignment in data.values():
        counter += 1
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
        conserved_regions += len(conserved_pos)
        matrix = update_matrix(matrix, conserved_pos, alignment_HCs, common_HCs)
    matrix = symmetrize(matrix)
    print('done')
    # save results
    np.savetxt(output_file, matrix, delimiter=' ', fmt='%d')
    return conserved_regions