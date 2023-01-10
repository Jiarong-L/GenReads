import configparser
import math
import pandas as pd
from GenReadsfunc.Genome_FileHandlers import Scan_GenomeSize


def load_param_cfg(param_cfg_dir):
    param_cfg = configparser.ConfigParser()
    param_cfg.read(param_cfg_dir)
    return param_cfg

def load_file_cfg(param_cfg):
    file_cfg_dir = param_cfg.get('inputs','file_cfg')
    merged_read_len = int(param_cfg.get('params','merged_read_len'))
    file_cfg = pd.read_csv(file_cfg_dir,sep='\t',index_col=0)
    assert file_cfg.index.is_unique
    assert set(file_cfg.columns) == set(['Amount', 'RawGenome'])
    file_cfg['merged_bases'] = file_cfg.apply(lambda x: convert_Size_Symbol(x[0],x[1]), axis = 1)
    file_cfg['merged_reads'] = file_cfg['merged_bases'].apply(lambda x: math.ceil(x/merged_read_len))
    return file_cfg.T.to_dict()



def convert_Size_Symbol(input_str,file_dir):
    amount_str = input_str.strip()
    size_unit = amount_str[-1].upper()
    if size_unit == 'X':
        total_bases = int(amount_str[:-1]) * Scan_GenomeSize(file_dir)
    elif size_unit == 'B':
        total_bases = int(amount_str[:-1])
    elif size_unit == 'K':
        total_bases = int(amount_str[:-1]) * 1024
    elif size_unit == 'M':
        total_bases = int(amount_str[:-1]) * 1024 * 1024
    elif size_unit == 'G':
        total_bases = int(amount_str[:-1]) * 1024 * 1024 *1024
    elif size_unit == 'T':
        total_bases = int(amount_str[:-1]) * 1024 * 1024 *1024 *1024
    elif size_unit in '0123456789':
        total_bases = int(amount_str)
    else:
        print("Size Symbol unsupported!")
        exit(0)
    return total_bases


