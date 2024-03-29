import numpy as np
transtab = str.maketrans('ATCGNSWRYMKHDBV','TAGCNSWYRKMDHVB')


def t_r(seq):
    return ''.join(reversed(seq.translate(transtab)))



def error_Base(temp_base):
    c = ['A','T','C','G','N']
    c.remove(temp_base)
    return np.random.choice(c,1)[0]

def Introduce_error(fragment,error_rate):
    if error_rate == 0:
        return fragment
    fragment = list(fragment)
    if_error_lst = np.random.choice([0,1],len(fragment),p=[1-error_rate,error_rate])
    error_idx = np.argwhere(if_error_lst == 1).squeeze()
    try:
        for i in error_idx:
            fragment[i] = error_Base(fragment[i])
    except:
        pass
    return ''.join(fragment)

def Split_R1R2(fragment,read_len):
    r1 = fragment[:read_len]
    r2_t_r = fragment[-read_len:]
    r2 = t_r(r2_t_r)
    return r1,r2


def Gen_Quality_line(r):
    return 'F'*len(r)



def Seq2Read(fragment,read_len,error_rate):
    r1,r2 = Split_R1R2(fragment,read_len)
    r1 = Introduce_error(r1,error_rate)
    r1_q = Gen_Quality_line(r1)
    r2 = Introduce_error(r2,error_rate)
    r2_q = Gen_Quality_line(r2)
    return r1,r1_q,r2,r2_q

def SaveSeq2Read(fragment,read_len,error_rate,header,fa_w,r1_w,r2_w):   # TODO: check if frag is valid & r1r2<frag
    r1,r1_q,r2,r2_q = Seq2Read(fragment,read_len,error_rate)
    fa_w.write_fa(header,fragment)
    r1_w.write_fq(header,r1,r1_q)
    r2_w.write_fq(header,r2,r2_q)

