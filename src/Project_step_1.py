import os
from time import time
import re
os.chdir("C:/Users/Alexis Trang/Documents/Cours_UPMC_M2/GENOM")


#####EXEMPLE
phrase="123 A:45-67"
phrase1="116748\tsp\ta.1.1.1\t-\tBacillus subtilis [TaxId: 1423]\n"
#############
class Domain:
    def __init__(self,ID,SCOPe,domain,PDB=None,chain=None,location=None,Pfam_accession=None):
        self.ID=ID
        self.SCOPe=SCOPe
        self.domain=domain
        self.PDB=PDB
        self.chain=chain
        self.location=location
        self.Pfam_accession=Pfam_accession
############################

soluble_domains=[]  #List of domains
PDBs=[]   #List of PDB codes of soluble domains

scope_file=open("dir.des.scope.2.08-stable.txt",'r')
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
            soluble_domains.append(Domain(line[0],line[2],line[3],PDB,chain,location))
            PDBs.append(PDB)
            soluble_domains_file.write(line[0]+"\t"+line[2]+"\t"+line[3]+"\t"+PDB+chain+str(location)+"\n")
    line=scope_file.readline()

scope_file.close()
soluble_domains_file.close()







