import os
from time import time

#PLEASE RUN Project_step_1.py before running that
os.chdir("C:/Users/Alexis Trang/Documents/Cours_UPMC_M2/GENOM")
pfam_mapping=open("pdb_pfam_mapping.txt",'r')

t1=time()
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
print(time()-t1)


############
"""
while line!="":
    line=line.split()
    if line[0] in PDBs:
        indices_PDBs=[ind for ind in range(len(PDBs)) if PDBs[ind]==line[0]]
        #print(indices_PDBs)
        #print("###################\n")
        for i in indices_PDBs:
            #print("sol:",soluble_domains[i].chain,"pfam:",line[1])

            if soluble_domains[i].chain==line[1]:
                #print("MATCH")
                print("sol:",soluble_domains[i].location,"pfam:","[\'"+line[2]+"-"+line[3]+"\']")
                if soluble_domains[i].location==line[2]+"-"+line[3]:
                    print("YUP")
                    soluble_domains[i].Pfam_accession=line[4]
                    PDBs[i]=-1
            #print(soluble_domains[i].PDB)
        #print("###################\n*\n*\n*\n*\n*\n*")

        #print(line[0],line[4])
    line=pfam_mapping.readline()
"""