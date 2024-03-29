

class Primers(object):
    def __init__(self):
        self.primer_dict = {}
        self.primer_dict['8F'] = 'AGAGTTTGATYMTGGCTCAG'
        self.primer_dict['27F'] = 'AGAGTTTGATCMTGGCTCAG'
        self.primer_dict['338F'] = 'ACTCCTACGGGAGGCAGCAG'
        self.primer_dict["341F"] = 'CCTACGGGNGGCWGCAG' 
        self.primer_dict["515F"] = 'GTGCCAGCMGCCGCGGTAA'
        self.primer_dict["515F_modified"] = 'GTGYCAGCMGCCGCGGTAA'
        self.primer_dict['518R'] = 'ATTACCGCGGCTGCTGG'
        self.primer_dict['806R'] = 'GGACTACHVGGGTWTCTAAT'
        self.primer_dict['806R_modified'] = 'GGACTACNVGGGTWTCTAAT'
        self.primer_dict['805R'] = 'GACTACHVGGGTATCTAATCC'
        self.primer_dict['909R'] = 'CCGTCAATTCMTTTRAGT'
        self.primer_dict['926R'] = 'CCGYCAATTYMTTTRAGTTT'
        self.primer_dict['1492R'] = 'AGAGTTTGATCMTGGCTCAG'
    def Get_Primer(self,primer):# http://journals.im.ac.cn/html/actamicrocn/2021/5/20210502.htm
        if primer not in self.primer_dict:    
            if self.Check_Valid(primer):
                self.primer_dict[primer] = primer
        return self.primer_dict[primer]
    def Check_Valid(self,primer): # Todo: check if primer is a valid seq
        return True


    