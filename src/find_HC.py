def binarization(seq):
    # list of contiguous strong hydrophobic amino acids
    list_aa = ['V', 'I', 'L', 'F', 'M', 'Y', 'W']
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
        if seq_aa[i]=='P':
            boundaries.add(i)
            j=1
            while i+j<len(seq_bin):
                if seq_bin[i+j]=='1':
                    break
                boundaries.add(i+j)
                j+=1
    return boundaries # complexity = (n-3)+4*#boundaries

def get_HC(seq, boundaries):
    structs = []
    struct = ''
    for i in range(len(seq)):
        if i not in boundaries:
            struct += seq[i]
        elif len(struct)>1:
            structs.append(struct)
            struct = ''
    return structs

def get_Q_codes(structs):
    Q_codes = []
    for struct in structs:
        Q_code =''
        for i in range(len(struct)-1):
            if struct[i:i+2]=='11':
                Q_code+='V'
            elif struct[i:i+3]=='101':
                Q_code+='M'
            elif struct[i:i+4]=='1001':
                Q_code+='U'
            elif struct[i:i+5]=='10001':
                Q_code+='D'
        if Q_code!='':
            Q_codes.append(Q_code)
    return Q_codes

def binary_coding(seq):
    # remove "." and "-"
    seq = seq.replace('.', '')
    seq = seq.replace('-', '')
    # transform into binary sequence
    seq_bin = binarization(seq)
    # find boundaries where there are not particular structure (helix or sheet)
    boundaries = find_intervals(seq_bin, seq)
    # get potential structures
    structs = get_HC(seq_bin, boundaries)
    # get Q-codes
    Q_codes = get_Q_codes(structs)
    
    return structs, Q_codes

def read_data(file):
    f = open(file, 'r')
    data = []
    data_i = []
    for line in f.readlines():
        if line[0]=='/':
            data.append(data_i)
            data_i=[]
        elif line[0]!='#':
            data_i.append(line.split()[1])
    f.close()
    return data

def count_HC(structs, dico):
    for struct in structs:
        for i in range(len(struct)-1):
            if struct[i]=='1':
                for j in range(i+1, len(struct)):
                    if struct[j]=='1':
                        if int(struct[i:j+1], 2) in dico:
                            dico[int(struct[i:j+1], 2)]+=1
                        else:
                            dico[int(struct[i:j+1], 2)] = 1
    return dico

def get_analyse(file):
    data = read_data(file)
    all_HC = {}
    all_Q_codes = []
    for i in range(len(data)):
        Q_codes = []
        for j in range(len(data[i])):
            structs, Q_code = binary_coding(data[i][j])
            all_HC = count_HC(structs, all_HC)
            Q_codes.append(Q_code)
        all_Q_codes.append(Q_codes)
    return all_Q_codes, all_HC