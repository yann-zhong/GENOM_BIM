import os
import numpy as np
from itertools import *
import re
from time import time

#############  Creation of the "Domain" class to easily store Domains and their informations
class Domain:
    def __init__(self,ID,SCOPe,domain,PDB=None,chain=None,location=None,Pfam_accession=None):
        self.ID=ID
        self.SCOPe=SCOPe
        self.domain=domain
        self.PDB=PDB
        self.chain=chain
        self.location=location
        self.Pfam_accession=Pfam_accession

############ Step 1: Soluble domains extraction

#SDE is a function that takes as input the path to a SCOPe file (str) and in output writes a new file that only contains soluble domains.

#The output file is called "soluble_domains.txt" and is created in the current directory.
#The output file is formated as such: ID	SCOPe	domain	PDB	chain	location
#The is also 2 variable outputs to this fuction:    SDE(path...) [0]: soluble_domains = every soluble domains under the class Domain created earlier in a list
#                                                   SDE(path...) [1]: PDBs =  the list of every PDB code for each domain in the soluble_domains list

#While running this pipeline, you should store those variables respectively under the names "soluble_domains" and "PDBs"
#For exemple: soluble_domains,PDBs=SDE("dir.des.scope.2.08-stable.txt")

def SDE(path_to_scope_file):
    soluble_domains=[]  #List of domains
    PDBs=[]   #List of PDB codes of soluble domains

    scope_file=open(path_to_scope_file,'r')
    soluble_domains_file=open("soluble_domains.txt",'w')

    for i in range(5):
        line=scope_file.readline()


    while(line!=""):
        line=line.split("\t")
        line[-1]=line[-1][:-1]
        if line[2][0] in ['a','b','c','d']:
            if re.findall(r" [A-Z]:",line[4])==[]:
                soluble_domains.append(Domain(line[0],line[2],line[3],line[4]))
                PDBs.append(line[4])
                soluble_domains_file.write(line[0]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\n")
            else:
                PDB=re.findall(r".+ ",line[4])[0][:-1]
                chain=re.findall(r"[A-Z]:",line[4])[0][0]
                location=re.findall(r"[0-9]+-[0-9]+",line[4])
                if location==[]:
                    location="full"
                else:
                    location=location[0]
                soluble_domains.append(Domain(line[0],line[2],line[3],PDB,chain,location))
                PDBs.append(PDB)
                soluble_domains_file.write(line[0]+"\t"+line[2]+"\t"+line[3]+"\t"+PDB+"\t"+chain+"\t"+str(location)+"\n")
        line=scope_file.readline()

    scope_file.close()
    soluble_domains_file.close()
    print("Soluble domains extracted in \"soluble_domains.txt\"")
    return soluble_domains,PDBs


############ Step 2: Pfam accession codes extraction

#PACE is a fuction that takes as input the path to a PDB/PFAM mapping file (str), a list of domains of interest (list of Domains), and the list of the domains' PDBs that should be sorted in the same order (list of str)

#The output file is called "soluble_domains_pfam_accession.txt" and is created in the current directory, it contains the same thing as the output SDE file but also includes the Pfam accession code for each domains.
#The output file is formated as such: ID	SCOPe	domain	PDB	chain	location Pfam_accession

def PACE(path_to_pdb_pfam_mapping_file,soluble_domains,PDBs):
    pfam_mapping=open(path_to_pdb_pfam_mapping_file,'r')
    Pfam_accessions=[]

    for i in range(2):
        line=pfam_mapping.readline()

    columns=line.split()

    line=pfam_mapping.readline()


    while line!="":

        soluble_domains_with_pfam_file=open("soluble_domains_pfam_accession.txt","a")
        line=line.split()
        if line[0] in PDBs:
            #print(line)
            indices_PDBs=[ind for ind in range(len(PDBs)) if PDBs[ind]==line[0]]
            #print(indices_PDBs)
            #print("###################\n")
            for i in indices_PDBs:
                #print("sol:",soluble_domains[i].chain,"pfam:",line[1])

                if soluble_domains[i].chain==line[1]:
                    if soluble_domains[i].location=='full':
                        soluble_domains_with_pfam_file.write(soluble_domains[i].ID+"\t"+soluble_domains[i].SCOPe+"\t"+soluble_domains[i].domain+"\t"+soluble_domains[i].PDB+"\t"+soluble_domains[i].chain+"\t"+line[2]+"-"+line[3]+"\t"+line[4]+"\n")
                        PDBs[i]=-1

                    else:
                        start=re.search(r"[\d]+-",soluble_domains[i].location).group()[:-1]
                        end=re.search(r"-[\d]+",soluble_domains[i].location).group()[1:]
                        #print(start,end)
                        if int(start)<=int(line[2])<int(line[3])<=int(end):
                            soluble_domains_with_pfam_file.write(soluble_domains[i].ID+"\t"+soluble_domains[i].SCOPe+"\t"+soluble_domains[i].domain+"\t"+soluble_domains[i].PDB+"\t"+soluble_domains[i].chain+"\t"+line[2]+"-"+line[3]+"\t"+line[4]+"\n")
                            PDBs[i]=-1
                    #print("MATCH")
            #print([ind for ind in range(len(PDBs)) if PDBs[ind]==line[0]])
        soluble_domains_with_pfam_file.close()
                #print(soluble_domains[i].PDB)
            #print("###################\n*\n*\n*\n*\n*\n*")

            #print(line[0],line[4])
        line=pfam_mapping.readline()

    soluble_domains_with_pfam_file.close()
    pfam_mapping.close()
    return "Soluble domains + Pfam accession codes extracted in \"soluble_domains_pfam_accession.txt\""


############ ADD THE REST OF THE PIPELINE???

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

############ Last step: Matrix construction

#MC is a function taht takes as input the path to a file that contains the names of all the selected HCs (str), and the path to the countinf matrix of theses chosen HCs (str)

#The output file is called "Matrix_swap.txt" and is created in the current directory, it contains the subtitution matrix of the given HCs.

def MC(path_all_HCs_names,path_count_matrix):

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
    return "Matrix created in \"Matrix_swap.txt\""





