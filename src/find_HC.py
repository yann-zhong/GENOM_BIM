import numpy as np

# Transform into binary sequence where the "1" represent the contiguous strong hydrophobic amino acids
# input: sequence of AA
# output: sequence in binary
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


# Find domains where there are no particular structures (helix or sheet)
# input: sequence in binary and sequence of AA
# output: set of positions not in particular structure
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
    return boundaries

# Get potential hydrophobic clusters
# input: sequence in binary and set of positions without particular structure
# output: list of HC limits and list of HC patterns
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

# Put back "." and "-" in the sequence
# input: list of HC limits, list of ponctuation positions
# output: list of adjusted HC limits
def put_back_ponctuation(positions, list_ponctuations):
    for ponct in list_ponctuations:
        for i in range(len(positions)):
            if positions[i]>=ponct:
                positions[i]+=1
    return positions

# Create a dictionary with positions of the HC and the HC
# input: list of HC limits and list of HC patterns
# output: dictionary of positions tuples (start, end) associated to their HC patterns
def make_dict(positions, HCs):
    final_dict = {}
    for i in range(int(len(positions)/2)):
        final_dict[(positions[i*2], positions[i*2+1])] = HCs[i]
    return final_dict

# Find the different HC patterns in a sequence and give their associated positions
# input: sequence of AA
# output: dictionary of positions tuples (start, end) associated to their HC patterns
def binary_coding(seq):
    # remove "." and "-" from the sequence but save their position
    list_ponctuations = [ind for ind,nucl in enumerate(seq) if nucl in ['.', '-']]
    seq = seq.replace('.', '').replace('-', '')
    seq_bin = binarization(seq)
    boundaries = find_intervals(seq_bin, seq)
    positions, HCs = get_HC(seq_bin, boundaries)
    positions = put_back_ponctuation(positions, list_ponctuations)
    return make_dict(positions, HCs)

# Get the Pfam accession numbers
# input: file containing the soluble domains of known 3D structures + Pfam accession on last column
# output: list of Pfam accession numbers
def get_PF(file):
    f = open(file, 'r')
    list_PF = set()
    for line in f.readlines():
        list_PF.add(line.split('\t')[-1].strip())
    f.close()
    return list_PF

# Read the alignments of soluble domains in a database
# input: file containing alignments and file containing soluble domains
# output: dictionary associating superfamily to their alignment
def read_data(file_alignments, file_PF):
    list_PF = get_PF(file_PF)
    f = open(file_alignments, 'r', encoding = 'ISO-8859-1')
    data = {}
    data_i = []
    count=1
    for line in f.readlines():
        if line[0]=='/':
            if superfamily in list_PF:
                data[superfamily]=data_i
                data_i=[]
            count=1
        elif line[0]!='#' and superfamily in list_PF:
            data_i.append(line.split()[1])
        elif count==3:
            superfamily = line.split()[2].split('.')[0]
            count+=1
        elif count<5:
            count+=1
    f.close()
    return data

# Count how many each HC patterns are present in the soluble domains and update the count
# input: dictionary of positions tuples (start, end) associated to their HC patterns and dictionary associating each HC pattern to its number of apparitions
# output: dictionary associating each HC pattern to its number of apparitions
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

# Get the Q-code of a HC pattern
# input: integer of a HC
# output: Q-code string
def get_Q_codes(hc):
    hc=bin(hc)
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

# Keep only the most present HC patterns
# input: dictionary associating HC patterns to its number of appartitions and the number of HC patterns kept
# output: adjusted dictionary associating HC patterns to its number of apparitions
def filtred_HCs(HCs, n):
    size = len(HCs)
    lim = int(np.quantile(list(HCs.values()), 1-n/size))
    HCs_filtred = {}
    for hc,count in HCs.items():
        if count>=lim:
            HCs_filtred[hc]=count
    return HCs_filtred

# Find most present HC patterns into soluble domain alignments and give potentially their Q-code
# input: file containing alignments, file containing soluble domains, number of HC patters kept and boolean which indicates if we need the Q-code
# output: dictionary associating HC patterns to his number of apparitions
def get_analyse(file_alignements, file_soluble_domains, nb_HC=500, Q_code=False):
    print('Preprocessing...')
    data = read_data(file_alignements, file_soluble_domains)
    print('done')
    print('Searching HCs...')
    all_HC = {}
    for _,alignment in data.items():
        for seq in alignment:
            HCs = binary_coding(seq)
            all_HC = count_HC(HCs, all_HC)
    print('done')
    print('Filtring...')
    all_HC = filtred_HCs(all_HC, nb_HC)
    all_HC = dict(sorted(all_HC.items(), key=lambda x:x[1], reverse=True))
    print('done')
    if Q_code:
        HC_with_Q_code = {}
        for hc,effectif in all_HC.items():
            HC_with_Q_code[(hc, get_Q_codes(hc))] = effectif
        return HC_with_Q_code
    return all_HC