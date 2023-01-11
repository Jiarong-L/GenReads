
from GenReads import *




cfg_file = 'cfg.ini'
MainCall(cfg_file)









# file_cfg.to_csv(os.path.join(saving_folder,'Gen_record(Random).xls'),sep='\t')





# ##############Tests for file handler

# file_dir = 'raw/fake.fna.gz'
# fh  = GenomeFastaReadHandler(file_dir)
# while fh.hasNextContig():
#     ContigHeader,ContigSeq = fh.nextContig()
#     print(ContigHeader)
#     print(ContigSeq)


# ts = Scan_GenomeSize(file_dir)
# print(ts)



# linkc = Link_Contigs(file_dir,500)
# print(linkc)
