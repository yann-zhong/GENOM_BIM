import numpy as np
from math import *

# Replace certain amino acids with sequences
# input: string sequence
# output: string sequence
def put_symbols(seq):
    dict_symbols = {
        'P':'★',
        'G':'⯁',
        'T':'☐',
        'S':'▣'
    }
    for aa in dict_symbols:
        seq = seq.replace(aa, dict_symbols[aa])
    return seq

# display the sequence in double helix form
# input: string sequence
# output: numpy matrix
def display_matrix(seq):
    n = len(seq)
    double_helix = np.empty((8, ceil(n/4)+1), dtype = "object")
    for i in range(n):
        line = i%4
        column = int(i/4)
        double_helix[line, 1+column] = seq[i]
        double_helix[line+4, column] = seq[i]
    double_helix[double_helix==None] = '-'
    print(double_helix)

# show the horizontal clusters
# input: string sequence
# output: numpy matrix
def horizontal_cluster(seq):
    # remove "." and "-"
    seq = (seq.replace('.', '')).replace('-', '')
    # transform into another form with special symbols denoting amino acids with particular structural behaviours
    seq = put_symbols(seq)
    # vizualize horizontal clusters
    display_matrix(seq)