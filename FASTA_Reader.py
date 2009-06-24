'''
Created on Jun 22, 2009

@author: Administrator
'''

class FASTA_Reader(object):
    '''
    class used to Extract sequences from FASTA db files
    '''
    File=None


    def __init__(selfparams,filename):
        '''
        Constructor
        '''
        selfparams.File=open(filename)

        
    def GetNextSequence(self):
        seq=''
        while True:
            line= self.File.readline()
            if line.startswith('>'):
                continue
            while not line.startswith('>') and len(line)>0:
                seq+=line.strip()
                line= self.File.readline()
            else:
                break
            if len(line)==0: break
        return seq
            
