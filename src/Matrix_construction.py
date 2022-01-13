import os
import numpy as np
from itertools import *
import re

############ Last step: Matrix construction

#MC is a function taht takes as input the path to a file that contains the names of all the selected HCs (path_all_HCs_names=str),  the path to the countinf matrix of theses chosen HCs (path_count_matrix=str), a scaling factor (lamda=float, if not specified, lambda=0.347), and a the path to a folder where you want your results saved (out_folder=str).

#The output file is called "Matrix_swap.txt" and is created in the current directory, it contains the subtitution matrix of the given HCs.

def MC(path_all_HCs_names,path_count_matrix,lamda=0.347,out_folder="."):

    names=np.loadtxt(path_all_HCs_names)#LOADING NAMES OF THE HCs (it's actually just the q_code for each HCs)
    n=len(names)  #NUMBER OF HCs

    matrice = np.loadtxt(path_count_matrix) #MATRICE LOADED

    dic_q={}  #INITIALIZING THE DICTIONARY WITH EVERY "qi" = MODEL NUL

    matrice_pij=matrice/np.sum(matrice)
    n = len(matrice)
    for i in range(n):
        if np.sum(matrice_pij[i])==0:
            matrice_pij[i][0]=1
        dic_q[int(names[i])]=(np.sum(matrice_pij[i])/2)+(matrice_pij[i,i])/2
    matrice_triangulaire = np.zeros(matrice.shape) #TRIANGLE

    for i in range(n):
        for j in range(n-(i+1)):
            matrice_triangulaire[i,j+(i+1)]=round(np.log(matrice_pij[i,j+(i+1)]/(2*dic_q[names[i]]*dic_q[names[j+(i+1)]]))/lamda)

    for i in range(n):
        matrice_triangulaire[i,i]=round(np.log(matrice_pij[i,i]/(dic_q[names[i]]**2))/lamda)

    matrice_file=open(out_folder+"\\Matrix_swap.txt","w")
    for i in matrice_triangulaire:
        for j in i:
            matrice_file.write(str(j)+" ")
        matrice_file.write("\n")
    return "Matrix created in \"Matrix_swap.txt\""



