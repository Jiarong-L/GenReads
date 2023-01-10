
from GenReadsfunc import *


# cfg_dir = 'raw/config_Random.xls'
# saving_folder='output'
# merged_read_length = 200

param_cfg_dir = 'param_cfg.ini'
param_cfg = load_param_cfg(param_cfg_dir)






file_cfg = load_file_cfg(param_cfg)
print(file_cfg)
# file_cfg.to_csv(os.path.join(saving_folder,'Gen_record(Random).xls'),sep='\t')





# ##############Tests for file handler

# file_dir = 'raw/fake.fna.gz'
# fh  = GenomeFastaReadHandler(file_dir)
# while fh.hasNextContig():
#     ContigHeader,ContigStr = fh.nextContig()
#     print(ContigHeader)
#     print(ContigStr)


# ts = Scan_GenomeSize(file_dir)
# print(ts)



# linkc = Link_Contigs(file_dir,500)
# print(linkc)
