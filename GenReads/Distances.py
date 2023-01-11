from GenReads.DegenerateDict import *
import numpy as np


degenerate_dict = Gen_degenerate_dict()


def blast():  ## TODO
    pass





def bio_hamming(seq1,seq2,degenerate_dict,allow_unequal_tail = True):  
    max_seq_len = max(len(seq1),len(seq2))
    min_seq_len = min(len(seq1),len(seq2))
    if not allow_unequal_tail:
        assert max_seq_len == min_seq_len
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


def Primer_Dist_scanner(ContigSeq,Primer,degenerate_dict): #,fill_marker,treashhold = 2
    contig_len = len(ContigSeq)
    Primer_len = len(Primer)
    dist_arr = np.array([])
    for i in range(contig_len):
        distance = bio_Levenshtein_dynamic(Primer,ContigSeq[i:i+Primer_len],degenerate_dict)# bio_hamming does not work!  use Levenshtein instead!
        dist_arr = np.append(dist_arr,distance)
    return dist_arr   # ,np.argwhere(dist_arr <= treashhold).squeeze(axis = 1) , min(dist_arr)

