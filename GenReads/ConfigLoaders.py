import configparser
import math
import pandas as pd
from GenReads.MainFunc import Scan_GenomeSize




def convert_Size_Symbol(input_str,ContigSizeList,insert_size):
    amount_str = input_str.strip()
    size_unit = amount_str[-1].upper()
    if size_unit == 'X':
        total_bases = int(amount_str[:-1]) * sum(ContigSizeList)
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
        total_bases = int(amount_str) * insert_size
    else:
        print("Size Symbol unsupported!")
        exit(0)
    return total_bases


def load_config_file(cfg_file):
    cfg_fh = configparser.ConfigParser()
    cfg_fh.read(cfg_file)
    cfg_dict = dict()
    for section in cfg_fh:
        for items in cfg_fh[section]:
            if cfg_fh[section][items] != '':
                cfg_dict[items] = cfg_fh[section][items]
    return cfg_dict


def load_several_input(cfg_dict):
    insert_size = int(cfg_dict['insert_size'])
    input_cfg = pd.read_csv(cfg_dict['input_cfg'],sep='\t',index_col=0)
    assert input_cfg.index.is_unique
    assert 'Amount' in set(input_cfg.columns)
    assert 'RawFasta' in set(input_cfg.columns)
    input_cfg['ContigSizeList'] = input_cfg.apply(lambda x: Scan_GenomeSize(x['RawFasta']), axis = 1)
    input_cfg['bases'] = input_cfg.apply(lambda x: convert_Size_Symbol(x['Amount'],x['ContigSizeList'],insert_size), axis = 1)
    input_cfg['frags'] = input_cfg['bases'].apply(lambda x: math.ceil(x/insert_size))
    return input_cfg.T.to_dict()


def load_single_input(cfg_dict):
    insert_size = int(cfg_dict['insert_size'])
    sample_dict = dict()
    sample_dict['Amount'] = cfg_dict['input_Amount'.lower()]
    sample_dict['RawFasta'] = cfg_dict['input_RawFasta'.lower()]
    sample_dict['ContigSizeList'] = Scan_GenomeSize(sample_dict['RawFasta'])
    sample_dict['bases'] = convert_Size_Symbol(sample_dict['Amount'],sample_dict['ContigSizeList'],insert_size)
    sample_dict['frags'] = math.ceil(sample_dict['bases']/insert_size)
    input_cfg_dict = {cfg_dict['input_SampleID'.lower()]:sample_dict}
    return input_cfg_dict


def load_input_info(cfg_dict):
    if 'input_cfg' in cfg_dict:   
        return load_several_input(cfg_dict)
    else:
        return load_single_input(cfg_dict)



