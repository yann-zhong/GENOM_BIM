import numpy as np
import os
from itertools import *

os.chdir("C:\\Users\Alexis Trang\\Documents\\Cours_UPMC_M2\\GENOM")

names=np.loadtxt("all_HCs.txt")#LOADING NAMES OF THE HCs (it's actually just the q_code for each HCs)
n=len(names)  #NUMBER OF HCs

matrice = np.loadtxt("matrix.txt") #MATRICE LOADED

dic_q={}  #INITIALIZING THE DICTIONARY WITH EVERY "qi" = MODEL NUL



matrice_pij=matrice/np.sum(matrice)

n = len(matrice)
for i in range(n):
    if np.sum(matrice_pij[i])==0:
        matrice_pij[i][0]=1
    dic_q[int(names[i])]=(np.sum(matrice_pij[i])/2)+(matrice_pij[i,i])/2


matrice_triangulaire = np.zeros(matrice.shape) #TRIANGLE

lamda=0.347



for i in range(n):
    for j in range(n-(i+1)):
        matrice_triangulaire[i,j+(i+1)]=round(np.log(matrice_pij[i,j+(i+1)]/(2*dic_q[names[i]]*dic_q[names[j+(i+1)]]))/lamda)

for i in range(n):
    matrice_triangulaire[i,i]=round(np.log(matrice_pij[i,i]/(dic_q[names[i]]**2))/lamda)

matrice_file=open("Matrix_swap.txt","w")
for i in matrice_triangulaire:
    for j in i:
        matrice_file.write(str(j)+" ")
    matrice_file.write("\n")








