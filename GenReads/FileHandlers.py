import gzip
import os



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
                        self.nextContigHeader = line[1:]
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
                    self.nextContigHeader = line[1:]
                    break 
                else:
                    ContigStr += line
        return ContigHeader,ContigStr.upper()




class WriterHandler(object):
    def __init__(self, file_dir):  
        self.file_dir = file_dir.strip()
        if os.path.exists(self.file_dir):      ## clean target file while init WriterHandler
            os.remove(self.file_dir)        
        if file_dir[-3:] == '.gz':        ## write to fa or gz file according to file suffix
            self.is_gzip = True
            self.fw = gzip.open(self.file_dir,'wb')
        else:
            self.is_gzip = False
            self.fw = open(self.file_dir,'w')

    def write_fq(self,header,seq,qualityline):
        fq_str = '@{}\n{}\n+\n{}\n'.format(header,seq,qualityline)
        if self.is_gzip:
            fq_str = fq_str.encode('utf-8')
        self.fw.write(fq_str)
        
    def write_fa(self,header,seq):
        fa_str = '>{}\n{}\n'.format(header,seq)
        if self.is_gzip:
            fa_str = fa_str.encode('utf-8')
        self.fw.write(fa_str)









