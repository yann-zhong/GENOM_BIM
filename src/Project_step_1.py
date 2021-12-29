import os
from time import time

os.chdir("C:/Users/Alexis Trang/Documents/Cours_UPMC_M2/GENOM")

class Domain:
    def __init__(self,ID,SCOPe,domain,chain):
        self.ID=ID
        self.SCOPe=SCOPe
        self.domain=domain
        self.chain=chain

############################

soluble_domains=[]

scope_file=open("dir.des.scope.2.08-stable.txt",'r')
for i in range(5):
    line=scope_file.readline()


while(line!=""):
    line=line.split("\t")
    line[-1]=line[-1][:-1]
    if line[2][0] in ['a','b','c','d']:
        soluble_domains.append(Domain(line[0],line[2],line[3],line[4]))
    line=scope_file.readline()



scope_file.close()







