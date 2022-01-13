import os
import numpy as np
from itertools import *
import re

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

#SDE is a function that takes as input the path to a SCOPe file (path_to_scope_file=str), and the path to a folder where you want your results saved (out_folder=str, if not specified, writes in the current directory).

#The output file is called "soluble_domains.txt" and only contains soluble domains.
#The output file is formated as such: ID	SCOPe	domain	PDB	chain	location
#There is also 2 variable outputs to this fuction:  SDE(path...) [0]: soluble_domains = every soluble domains under the class Domain created earlier in a list
#                                                   SDE(path...) [1]: PDBs =  the list of every PDB code for each domain in the soluble_domains list

#While running this pipeline, you should store those variables respectively under the names "soluble_domains" and "PDBs"
#For exemple: soluble_domains,PDBs=SDE("dir.des.scope.2.08-stable.txt")

def SDE(path_to_scope_file,out_folder="."):
    soluble_domains=[]  #List of domains
    PDBs=[]   #List of PDB codes of soluble domains

    scope_file=open(path_to_scope_file,'r')
    soluble_domains_file=open(out_folder+"\\soluble_domains.txt",'w')

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
    print("Soluble domains extracted in \""+out_folder+"/soluble_domains.txt\"")
    return soluble_domains,PDBs


############ Step 2: Pfam accession codes extraction

#PACE is a fuction that takes as input the path to a PDB/PFAM mapping file (path_to_pdb_pfam_mapping_file=str), a list of domains of interest (soluble_domains=list of Domains), the list of the domains' PDBs that should be sorted in the same order (PDBs=list of str),and the path to a folder where you want your results saved (out_folder=str, if not specified, writes in the current directory).

#The output file is called "soluble_domains_pfam_accession.txt", and it contains the same thing as the output SDE file but also includes the Pfam accession code for each domains.
#The output file is formated as such: ID	SCOPe	domain	PDB	chain	location Pfam_accession

def PACE(path_to_pdb_pfam_mapping_file,soluble_domains,PDBs,out_folder="."):
    pfam_mapping=open(path_to_pdb_pfam_mapping_file,'r')
    Pfam_accessions=[]

    for i in range(2):
        line=pfam_mapping.readline()

    columns=line.split()

    line=pfam_mapping.readline()


    while line!="":

        soluble_domains_with_pfam_file=open(out_folder+"\\soluble_domains_pfam_accession.txt","a")
        line=line.split()
        if line[0] in PDBs:
            indices_PDBs=[ind for ind in range(len(PDBs)) if PDBs[ind]==line[0]]
            for i in indices_PDBs:
                if soluble_domains[i].chain==line[1]:
                    if soluble_domains[i].location=='full':
                        soluble_domains_with_pfam_file.write(soluble_domains[i].ID+"\t"+soluble_domains[i].SCOPe+"\t"+soluble_domains[i].domain+"\t"+soluble_domains[i].PDB+"\t"+soluble_domains[i].chain+"\t"+line[2]+"-"+line[3]+"\t"+line[4]+"\n")
                        PDBs[i]=-1

                    else:
                        start=re.search(r"[\d]+-",soluble_domains[i].location).group()[:-1]
                        end=re.search(r"-[\d]+",soluble_domains[i].location).group()[1:]
                        if int(start)<=int(line[2])<int(line[3])<=int(end):
                            soluble_domains_with_pfam_file.write(soluble_domains[i].ID+"\t"+soluble_domains[i].SCOPe+"\t"+soluble_domains[i].domain+"\t"+soluble_domains[i].PDB+"\t"+soluble_domains[i].chain+"\t"+line[2]+"-"+line[3]+"\t"+line[4]+"\n")
                            PDBs[i]=-1
        soluble_domains_with_pfam_file.close()
        line=pfam_mapping.readline()

    soluble_domains_with_pfam_file.close()
    pfam_mapping.close()
    return "Soluble domains + Pfam accession codes extracted in \""+out_folder+"/soluble_domains_pfam_accession.txt\""



