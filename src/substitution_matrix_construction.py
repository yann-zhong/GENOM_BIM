import os
import numpy as np
from itertools import *
import re
import json
from src.count_matrix_construction import symmetrize

############ Last step: Matrix construction

#MC is a function that takes as input the names of all the selected HCs (names=list of str),  the path to the countinf matrix of theses chosen HCs 
# (path_count_matrix=str), a scaling factor (lamda=float, if not specified, lambda=0.347), 
# and a the path to a folder where you want your results saved (out_folder=str).

#The output file is called "Matrix_swap.txt" and is created in the current directory, it contains the subtitution matrix of the given HCs.

def MC(HC_counts, path_count_matrix, output_file, lamda=0.347):
    
    file = open(HC_counts, 'r')
    names = list(json.load(file).keys())
    file.close()
    n=len(names)  #NUMBER OF HCs
    for i in range(n):
        names[i]=int(names[i])

    matrix = np.loadtxt(path_count_matrix) #MATRICE LOADED

    dic_q={}  #INITIALIZING THE DICTIONARY WITH EVERY "qi" = MODEL NUL

    matrix_pij=matrix/np.sum(matrix)
    
    n = len(matrix)
    for i in range(n):
        dic_q[names[i]]=(np.sum(matrix_pij[i])/2)+(matrix_pij[i,i])/2

    matrix_triangulaire = np.zeros(matrix.shape) #TRIANGLE

    for i in range(n-1):
        for j in range(i+1, n):
            matrix_triangulaire[i,j]=round(np.log(matrix_pij[i,j]/(2*dic_q[names[i]]*dic_q[names[j]]))/lamda)

    for i in range(n):
        matrix_triangulaire[i,i]=round(np.log(matrix_pij[i,i]/(dic_q[names[i]]**2))/lamda)

    matrix_sym = symmetrize(matrix_triangulaire)
    
    matrix_file=open(output_file, "w")
    for i in matrix_sym:
        for j in i:
            matrix_file.write(str(j)+" ")
        matrix_file.write("\n")
    return "Matrix created in"+output_file