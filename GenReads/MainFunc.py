from GenReads.FileHandlers import *
from GenReads.Distances import *
from GenReads.Primers import *
from GenReads.ReadProcessors import *
import numpy as np




def Scan_GenomeSize(file_dir):
    size_counter = 0
    f = GenomeFastaReadHandler(file_dir)
    while f.hasNextContig():
        _,ContigStr = f.nextContig()
        size_counter += len(ContigStr)
    return size_counter



def Link_Contigs(file_dir,gapSize,gapfiller='N'): ## gapSize = insert_size-1     ## oldname: load_Genome
    gap_seq = gapfiller*gapSize
    linked_Genomes = ''
    f = GenomeFastaReadHandler(file_dir)
    while f.hasNextContig():
        _,ContigStr = f.nextContig()
        if linked_Genomes != '':
            linked_Genomes += gap_seq
        linked_Genomes += ContigStr
    return linked_Genomes



def extract_amplicon(ContigSeq,Primer_F,Primer_R_tr,degenerate_dict,insert_size,insert_size_max,insert_size_min,size_punish=0.02,F_treash = 3,R_treash = 3):
    F_dist_arr = Primer_Dist_scanner(ContigSeq,Primer_F,degenerate_dict)
    R_dist_arr = Primer_Dist_scanner(ContigSeq,Primer_R_tr,degenerate_dict)
    res_lst = []
    punish_arr = np.array([])
    for F in np.argwhere(F_dist_arr <= F_treash).squeeze(axis = 1):
        for R in np.argwhere(R_dist_arr <= R_treash).squeeze(axis = 1):
            res_seq = ContigSeq[F:R+len(Primer_R_tr)+1]
            if (len(res_seq) <= insert_size_max):
                if (len(res_seq) >= insert_size_min):
                    res_lst.append(res_seq)
                    punish = abs(R+len(Primer_R_tr)-F+1 - insert_size)*size_punish + F_dist_arr[F] + R_dist_arr[R]
                    punish_arr = np.append(punish_arr, punish)
    return punish_arr,res_lst


def AnchorModeGeneration():  ## TODO: anchor pcr   single primer?
    pass



def PrimerModeGeneration(ContigSeq,cfg_dict,sample_dict,prefix,fa_w,r1_w,r2_w):
    Primer_F = Primers().Get_Primer(cfg_dict['forward'])
    Primer_R = Primers().Get_Primer(cfg_dict['reverse'])
    Primer_R_tr = t_r(Primer_R)
    insert_size = int(cfg_dict['insert_size'])
    insert_size_max = int(cfg_dict['insert_size_max'])
    insert_size_min = int(cfg_dict['insert_size_min'])
    read_len = int(cfg_dict['read_len'])
    error_rate = float(cfg_dict['error_rate'])

    punish_arr, res_lst = extract_amplicon(ContigSeq,Primer_F,Primer_R_tr,degenerate_dict,
                                            insert_size,insert_size_max,insert_size_min,
                                            size_punish=0.02,F_treash = 3,R_treash = 3)
    # TODO: prob_arr = 1 / punish_arr
    if len(res_lst)>0:
        counter = 1
        for idx in np.argwhere(punish_arr <= min(punish_arr)).squeeze(axis=1):  # Select the most possible result
            header = '{}__{}'.format(prefix,counter)
            counter += 1
            best_frag = res_lst[idx][len(Primer_F):-len(Primer_R_tr)]           # drop primer seq TODO:fix it
            SaveSeq2Read(best_frag,read_len,error_rate,header,fa_w,r1_w,r2_w)



def RandomModeGeneration(ContigSeq,cfg_dict,sample_dict,prefix,fa_w,r1_w,r2_w):
    insert_size = int(cfg_dict['insert_size'])
    total_reads = int(sample_dict['reads'])
    read_len = int(cfg_dict['read_len'])
    error_rate = float(cfg_dict['error_rate'])
    pos_lst = np.random.choice((len(ContigSeq)-insert_size),total_reads) #replace=True, can be selected multiple times
    for i in range(len(pos_lst)):  
        header = '{}__{}'.format(prefix,i)
        random_pos =  pos_lst[i]
        random_frag = ContigSeq[random_pos:random_pos+insert_size]
        assert len(random_frag) == insert_size
        SaveSeq2Read(random_frag,read_len,error_rate,header,fa_w,r1_w,r2_w)




def CallGenerators(ContigSeq,cfg_dict,sample_dict,prefix,fa_w,r1_w,r2_w):
    if 'forward' in cfg_dict:
        if 'reverse' in cfg_dict:
            return PrimerModeGeneration(ContigSeq,cfg_dict,sample_dict,prefix,fa_w,r1_w,r2_w)
        else:
            return AnchorModeGeneration()  ## TODO
    else:
        return RandomModeGeneration(ContigSeq,cfg_dict,sample_dict,prefix,fa_w,r1_w,r2_w)








