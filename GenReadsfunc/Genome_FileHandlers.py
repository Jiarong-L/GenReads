import gzip



class GenomeFastaReadHandler(object):
    def __init__(self, file_dir):  
        self.file_dir = file_dir.strip()
        self.is_FirstContig = True
        self.has_NextLine = True
        if file_dir[-3:] == '.gz':        ## load fa or gz file according to file suffix
            self.is_gzip = True
            self.fr = gzip.open(self.file_dir,'rb')
        else:
            self.is_gzip = False
            self.fr = open(self.file_dir,'r')

    def nextLine(self):                   ## fetch next line and update has_NextLine status
        if self.is_gzip:
            line = self.fr.readline().decode('utf-8')
        else:
            line = self.fr.readline()
        if not line:
            self.has_NextLine = False
            self.fr.close()              ## close filehandle
        return line.strip()

    def firstContigHeader(self):                ## loading first contig's header [call by hasNextContig()]
        if self.is_FirstContig:
            self.is_FirstContig = False
            while self.has_NextLine:
                line = self.nextLine()
                if (len(line)>0):
                    if line[0] == '>':
                        self.nextContigHeader = line
                        break

    def hasNextContig(self):
        self.firstContigHeader()
        if not self.nextContigHeader:
            return False
        return True

    def nextContig(self):                ## loading string of next contig and header of next*2 contig [if hasNextContig]
        if not self.hasNextContig():
            return None,None
        ContigStr = ''
        ContigHeader = self.nextContigHeader
        self.nextContigHeader = None
        while self.has_NextLine:
            line = self.nextLine()
            if (len(line)>0):
                if line[0] == '>':
                    self.nextContigHeader = line
                    break 
                else:
                    ContigStr += line
        return ContigHeader,ContigStr.upper()

    # def scanSize(self):
    #     size_counter = 0
    #     while self.has_NextLine:
    #         line = self.nextLine()
    #         if (len(line)>0):
    #             if line[0] == '>':
    #                 pass
    #             else:
    #                 size_counter += len(line)
    #     return size_counter




def Scan_GenomeSize(file_dir):
    size_counter = 0
    f = GenomeFastaReadHandler(file_dir)
    while f.hasNextContig():
        _,ContigStr = f.nextContig()
        size_counter += len(ContigStr)
    return size_counter



def Link_Contigs(file_dir,gapSize,gapfiller='N'): ## gapSize = merged_read_len-1     ## oldname: load_Genome
    gap_seq = gapfiller*gapSize
    linked_Genomes = ''
    f = GenomeFastaReadHandler(file_dir)
    while f.hasNextContig():
        _,ContigStr = f.nextContig()
        if linked_Genomes != '':
            linked_Genomes += gap_seq
        linked_Genomes += ContigStr
    return linked_Genomes











