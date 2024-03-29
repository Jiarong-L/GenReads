from GenReads.ConfigLoaders import *
from GenReads.FileHandlers import *
from GenReads.MainFunc import *



def MainCall(cfg_file):
    cfg_dict = load_config_file(cfg_file)
    input_cfg_dict = load_input_info(cfg_dict)
    link_contigs = cfg_dict['each_contig_as_a_genome'].upper() in ['False'.upper(),'F','No'.upper()]
    for sid in input_cfg_dict:
        fa_w = WriterHandler(os.path.join(cfg_dict['output_folder'],'{}.fa'.format(sid)))
        r1_w = WriterHandler(os.path.join(cfg_dict['output_folder'],'{}_R1.fq'.format(sid)))
        r2_w = WriterHandler(os.path.join(cfg_dict['output_folder'],'{}_R2.fq'.format(sid)))
        sample_dict = input_cfg_dict[sid]
        file_dir = sample_dict['RawFasta']
        if link_contigs:
            gapSize = int(cfg_dict['insert_size']) -1
            ContigSeq = Link_Contigs(file_dir,gapSize,gapfiller='N')
            CallGenerators(ContigSeq,cfg_dict,sample_dict,sid,'',fa_w,r1_w,r2_w)
        else:
            # sample_dict['1X_bases'] = 1  ## TODO   
            fh  = GenomeFastaReadHandler(file_dir)
            while fh.hasNextContig():
                ContigHeader,ContigSeq = fh.nextContig()
                CallGenerators(ContigSeq,cfg_dict,sample_dict,sid,ContigHeader,fa_w,r1_w,r2_w)

