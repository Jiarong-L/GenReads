import numpy as np
transtab = str.maketrans('ATCGNSWRYMKHDBV','TAGCNSWYRKMDHVB')

def Gen_degenerate_dict(fill_marker_key = ['-'],fill_marker_val = ['']):
    degenerate_dict = {
        'A':'A',
        'T':'T',
        'C':'C',
        'G':'G',
        'R':'AG',
        'Y':'CT',
        'M':'AC',
        'K':'GT',
        'S':'GC',
        'W':'AT',
        'H':'ATC',
        'B':'GTC',
        'V':'GAC',
        'D':'GAT',
        'N':'ATCG'
    }
    for key,val in zip(fill_marker_key,fill_marker_val):
        degenerate_dict[key] = val
    return degenerate_dict


def Gen_Primer(primer_name):# http://journals.im.ac.cn/html/actamicrocn/2021/5/20210502.htm
    primer_dict = {}
    primer_dict['8F'] = 'AGAGTTTGATYMTGGCTCAG'
    primer_dict['27F'] = 'AGAGTTTGATCMTGGCTCAG'
    primer_dict['338F'] = 'ACTCCTACGGGAGGCAGCAG'
    primer_dict["341F"] = 'CCTACGGGNGGCWGCAG' 
    primer_dict["515F"] = 'GTGCCAGCMGCCGCGGTAA'
    primer_dict["515F_modified"] = 'GTGYCAGCMGCCGCGGTAA'
    primer_dict['518R'] = 'ATTACCGCGGCTGCTGG'
    primer_dict['806R'] = 'GGACTACHVGGGTWTCTAAT'
    primer_dict['806R_modified'] = 'GGACTACNVGGGTWTCTAAT'
    primer_dict['805R'] = 'GACTACHVGGGTATCTAATCC'
    primer_dict['909R'] = 'CCGTCAATTCMTTTRAGT'
    primer_dict['926R'] = 'CCGYCAATTYMTTTRAGTTT'
    primer_dict['1492R'] = 'AGAGTTTGATCMTGGCTCAG'
    return primer_dict[primer_name]

def bio_hamming(seq1,seq2,degenerate_dict,allow_unequal_tail = True,ignore_case=True):  
    max_seq_len = max(len(seq1),len(seq2))
    min_seq_len = min(len(seq1),len(seq2))
    if not allow_unequal_tail:
        assert max_seq_len == min_seq_len
    if ignore_case:
        seq1 = seq1.upper()
        seq2 = seq2.upper()
    distance = 0
    for i in range(min_seq_len):
        if len(set(degenerate_dict[seq1[i]]).intersection(degenerate_dict[seq2[i]])) == 0:
            distance += 1
    distance += max_seq_len-min_seq_len
    return max_seq_len,distance


def bio_Levenshtein(seq1,seq2,degenerate_dict,insert_cost = 1,delet_cost = 1,mis_cost = 1): # too slow for long seq
    i = len(seq1)
    j = len(seq2)
    if min(i,j) == 0:
        return max(i,j)
    insertion_L = bio_Levenshtein(seq1,seq2[:-1],degenerate_dict,insert_cost = insert_cost,delet_cost = delet_cost,mis_cost = mis_cost) + insert_cost
    deletion_L = bio_Levenshtein(seq1[:-1],seq2,degenerate_dict,insert_cost = insert_cost,delet_cost = delet_cost,mis_cost = mis_cost) + delet_cost
    mismatch_L = bio_Levenshtein(seq1[:-1],seq2[:-1],degenerate_dict,insert_cost = insert_cost,delet_cost = delet_cost,mis_cost = mis_cost)
    if len(set(degenerate_dict[seq1[-1]]).intersection(degenerate_dict[seq2[-1]])) == 0:
        mismatch_L += mis_cost
    return min(insertion_L,deletion_L,mismatch_L)

def bio_Levenshtein_dynamic(seq1,seq2,degenerate_dict,insert_cost = 1,delet_cost = 1,mis_cost = 1):
    dim1 = len(seq1) +1 
    dim2 = len(seq2) +1
    dp = np.zeros((dim1,dim2))
    for i in range(dim1):
        dp[i][0] = i
    for j in range(dim2):
        dp[0][j] = j
    for i in range(1,dim1):
        for j in range(1,dim2):
            if_mis_cost = 0 if len(set(degenerate_dict[seq1[i-1]]).intersection(degenerate_dict[seq2[j-1]])) != 0 else 1
            mismatch_L = dp[i-1][j-1]+ if_mis_cost*mis_cost
            insertion_L = dp[i][j-1] + insert_cost
            deletion_L = dp[i-1][j] + delet_cost
            dp[i][j] = min( mismatch_L,insertion_L,deletion_L)
    return int(dp[i][j])


def scan_Primer(contig_seq,Primer,degenerate_dict): #,fill_marker,treashhold = 2
    contig_len = len(contig_seq)
    Primer_len = len(Primer)
    dist_arr = np.array([])
    for i in range(contig_len):
        distance = bio_Levenshtein_dynamic(Primer,contig_seq[i:i+Primer_len],degenerate_dict)# bio_hamming does not work!  use Levenshtein instead!
        dist_arr = np.append(dist_arr,distance)
    return dist_arr   # ,np.argwhere(dist_arr <= treashhold).squeeze(axis = 1) , min(dist_arr)


def extract_amplicon(contig_seq,Primer_F,Primer_R_tr,degenerate_dict,Amplicon_size,Amplicon_size_max,Amplicon_size_min,size_punish=0.02,F_treash = 3,R_treash = 3):
    F_dist_arr = scan_Primer(contig_seq,Primer_F,degenerate_dict)
    R_dist_arr = scan_Primer(contig_seq,Primer_R_tr,degenerate_dict)
    res_lst = []
    punish_arr = np.array([])
    for F in np.argwhere(F_dist_arr <= F_treash).squeeze(axis = 1):
        for R in np.argwhere(R_dist_arr <= R_treash).squeeze(axis = 1):
            res_seq = contig_seq[F:R+len(Primer_R_tr)+1]
            if (len(res_seq) <= Amplicon_size_max):
                if (len(res_seq) >= Amplicon_size_min):
                    res_lst.append(res_seq)
                    punish = abs(R+len(Primer_R_tr)-F+1 - Amplicon_size)*size_punish + F_dist_arr[F] + R_dist_arr[R]
                    punish_arr = np.append(punish_arr, punish)
    return punish_arr,res_lst



def formatting_contigs(file_dir,out_dir):
    with open(file_dir,'r') as fr,open(out_dir,'w') as fw:
        first_line = True
        while True:
            line = fr.readline()
            if not line:
                fw.write('\n')
                break
            if line[0] == '>':
                if first_line:
                    fw.write(line)
                    first_line = False
                else:
                    fw.write('\n')
                    fw.write(line)
            else:
                fw.write(line.strip())
    return True


def load_contigs(file_dir,contig_prefix):
    contig_dict = {}
    contig_counter = 0
    with open(file_dir,'r') as f:
        while True:
            line = f.readline().strip()
            if not line:
                break
            if line[0] == '>':
                contig_counter +=1
                contig_ID = '{}_{}_{}'.format(contig_prefix,contig_counter,line.strip('>'))
                contig_dict[contig_ID] = ''
            else:
                contig_dict[contig_ID] += line
    return contig_dict,contig_counter


# def load_single_contig(file_handle,contig_prefix,contig_counter):
#     if_next = True
#     contig_seq = ''
#     while True:
#         line = file_handle.readline().strip()
#         if not line:
#             if_next = False
#             return '',contig_seq,if_next
#         if line[0] == '>':
#             next_contig_ID = '{}_{}_{}'.format(contig_prefix,contig_counter,line.strip('>'))
#             return next_contig_ID,contig_seq,if_next
#         else:
#             contig_seq += line

# def load_contigs(file_dir,contig_prefix):  # for testing load_single_contig
#     contig_dict = {}
#     contig_counter = 1
#     with open(file_dir,'r') as f:
#         contig_ID,_,_ = load_single_contig(f,contig_prefix,contig_counter)
#         while True:
#             next_contig_ID,contig_seq,if_next = load_single_contig(f,contig_prefix,contig_counter)
#             if not if_next:
#                 break
#             else:
#                 contig_counter +=1 
#             contig_dict[contig_ID] = contig_seq
#             contig_ID = next_contig_ID
#     return contig_dict,contig_counter


def load_Genome(file_dir,merged_read_length):
    gap_seq = 'N'*(merged_read_length-1)
    with open(file_dir,'r') as f:
        linked_Genome = ''
        while True:
            line = f.readline().strip()
            if not line:
                break
            if line[0] == '>':
                linked_Genome += gap_seq
            else:
                linked_Genome += line
    return linked_Genome.upper()

def error_Base(temp_base):
    c = ['A','T','C','G','N']
    c.remove(temp_base)
    return np.random.choice(c,1)[0]

def Introduce_error(random_read,error_rate = 0.01):
    if error_rate == 0:
        return random_read
    random_read = list(random_read)
    if_error_lst = np.random.choice([0,1],len(random_read),p=[1-error_rate,error_rate])
    error_idx = np.argwhere(if_error_lst == 1).squeeze()
    try:
        for i in error_idx:
            random_read[i] = error_Base(random_read[i])
    except:
        pass
    return ''.join(random_read)

def Split_R1R2(random_read,single_read_length):
    r1 = random_read[:single_read_length]
    r2_t_r = random_read[-single_read_length:]
    r2 = ''.join(reversed(r2_t_r.translate(transtab)))
    return r1,r2


def write_fq(file_handle,header,seq,quality = 'F'):
    file_handle.write('@{}\n'.format(header))
    file_handle.write('{}\n'.format(seq))
    file_handle.write('+\n')
    file_handle.write('{}\n'.format(quality*len(seq)))
