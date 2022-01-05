import numpy as np

def binarization(seq):
    # list of contiguous strong hydrophobic amino acids
    list_aa = ['V', 'I', 'L', 'F', 'M', 'Y', 'W', 'v', 'i', 'l', 'f', 'm', 'y', 'w']
    output = ''
    for i in range(len(seq)):
        if seq[i] in list_aa:
            output += '1'
        else:
            output += '0'
    return output

def find_intervals(seq_bin, seq_aa):
    boundaries = set()
    for i in range(len(seq_bin)-3):
        if seq_bin[i:i+4]=='0000':
            for j in range(4):
                boundaries.add(i+j)
    for i in range(len(seq_aa)):
        if seq_aa[i] in ['P', 'p']:
            boundaries.add(i)
            j=1
            while i+j<len(seq_bin):
                if seq_bin[i+j]=='1':
                    break
                boundaries.add(i+j)
                j+=1
    return boundaries # complexity = (n-3)+4*#boundaries

def get_HC(seq, boundaries):
    positions = []
    HCs = []
    hc = ''
    for i in range(len(seq)):
        if i not in boundaries:
            if hc=='':
                pos1=i
            hc += seq[i]
        elif len(hc)>1:
            HCs.append(hc)
            positions += [pos1, i-1]
            hc = ''
        else:
            hc = ''
    if len(hc)>1:
        HCs.append(hc)
        positions += [pos1, len(seq)-1]
    return positions, HCs

def put_back_ponctuation(positions, list_ponctuations):
    for ponct in list_ponctuations:
        for i in range(len(positions)):
            if positions[i]>=ponct:
                positions[i]+=1
    return positions

def make_dict(positions, HCs):
    final_dict = {}
    for i in range(int(len(positions)/2)):
        final_dict[(positions[i*2], positions[i*2+1])] = HCs[i]
    return final_dict

def binary_coding(seq):
    # remove "." and "-" but conserve their position
    list_ponctuations = [ind for ind,nucl in enumerate(seq) if nucl in ['.', '-']]
    seq = seq.replace('.', '').replace('-', '')
    # transform into binary sequence
    seq_bin = binarization(seq)
    # find boundaries where there are not particular structure (helix or sheet)
    boundaries = find_intervals(seq_bin, seq)
    # get potential hydrophobic clusters
    positions, HCs = get_HC(seq_bin, boundaries)
    # put back "." and "-"
    positions = put_back_ponctuation(positions, list_ponctuations)
    # make a dictionnary with the positions and the HC
    return make_dict(positions, HCs)

def read_data(file):
    f = open(file, 'r', encoding = 'ISO-8859-1')
    data = {}
    data_i = []
    count=1
    for line in f.readlines():
        if line[0]=='/':
            data[superfamily]=data_i
            data_i=[]
            count=1
        elif line[0]!='#':
            data_i.append(line.split()[1])
        elif count==3:
            superfamily = line.split()[2].split('.')[0]
            count+=1
        elif count<5:
            count+=1
    f.close()
    return data

def get_PF(file):
    f = open(file)
    list_PF = set()
    for line in f.readlines():
        list_PF.add(line.split('\t')[-1].strip())
    f.close()
    return list_PF

def count_HC(HCs, dico):
    for hc in HCs.values():
        for i in range(len(hc)-1):
            if hc[i]=='1':
                for j in range(i+1, len(hc)):
                    if hc[j]=='1':
                        if int(hc[i:j+1], 2) in dico:
                            dico[int(hc[i:j+1], 2)]+=1
                        else:
                            dico[int(hc[i:j+1], 2)] = 1
    return dico

def get_Q_codes(hc):
    Q_code =''
    for i in range(len(hc)-1):
        if hc[i:i+2]=='11':
            Q_code+='V'
        elif hc[i:i+3]=='101':
            Q_code+='M'
        elif hc[i:i+4]=='1001':
            Q_code+='U'
        elif hc[i:i+5]=='10001':
            Q_code+='D'
    return Q_code

def filtred_HCs(HCs, n):
    size = len(HCs)
    lim = int(np.quantile(list(HCs.values()), 1-n/size))
    HCs_filtred = {}
    for hc,count in HCs.items():
        if count>=lim:
            HCs_filtred[hc]=count
    return HCs_filtred

def get_analyse(file_alignements, file_soluble_domains, nb_HC=500, Q_code=False):
    data = read_data(file_alignements)
    list_PF = get_PF(file_soluble_domains)
    list_PF = data.keys()
    all_HC = {}
    for superfamily,alignment in data.items():
        if superfamily in list_PF:
            for seq in alignment:
                HCs = binary_coding(seq)
                all_HC = count_HC(HCs, all_HC)
    all_HC = filtred_HCs(all_HC, nb_HC)
    all_HC = dict(sorted(all_HC.items(), key=lambda x:x[1], reverse=True))
    if Q_code:
        HC_with_Q_code = {}
        for hc in all_HC:
            HC_with_Q_code[hc] = get_Q_codes(hc)
        return HC_with_Q_code
    return all_HC