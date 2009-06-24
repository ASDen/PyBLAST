'''
Created on Jun 22, 2009

@author: Administrator
'''
import itertools
import time

class QueryProcessor(object):
    '''
    A class , used for pre-processing BLAST queries , so as to
    improve performance
    '''
    Prtn_Alpha = ['A' , 'R' , 'N' , 'D' , 'C' , 'E' , 'Q' , 'G' , 'H' , 'I' , 'L' , 'K' , 'M' , 'F' , 'P' , 'S' , 'T' , 'W' , 'Y' , 'V']
    
    Prtn_Alpha_dict={}
    Prtn_Alpha_length=0
    
    blosum62 = [
4 , -1 , -2 , -2 , 0 , -1 , -1 , 0 , -2 , -1 , -1 , -1 , -1 , -2 , -1 , 1 , 0 , -3 , -2 , 0
, -1 , 5 , 0 , -2 , -3 , 1 , 0 , -2 , 0 , -3 , -2 , 2 , -1 , -3 , -2 , -1 , -1 , -3 , -2 , -3
, -2 , 0 , 6 , 1 , -3 , 0 , 0 , 0 , 1 , -3 , -3 , 0 , -2 , -3 , -2 , 1 , 0 , -4 , -2 , -3
, -2 , -2 , 1 , 6 , -3 , 0 , 2 , -1 , -1 , -3 , -4 , -1 , -3 , -3 , -1 , 0 , -1 , -4 , -3 , -3
, 0 , -3 , -3 , -3 , 9 , -3 , -4 , -3 , -3 , -1 , -1 , -3 , -1 , -2 , -3 , -1 , -1 , -2 , -2 , -1
, -1 , 1 , 0 , 0 , -3 , 5 , 2 , -2 , 0 , -3 , -2 , 1 , 0 , -3 , -1 , 0 , -1 , -2 , -1 , -2
, -1 , 0 , 0 , 2 , -4 , 2 , 5 , -2 , 0 , -3 , -3 , 1 , -2 , -3 , -1 , 0 , -1 , -3 , -2 , -2
, 0 , -2 , 0 , -1 , -3 , -2 , -2 , 6 , -2 , -4 , -4 , -2 , -3 , -3 , -2 , 0 , -2 , -2 , -3 , -3
, -2 , 0 , 1 , -1 , -3 , 0 , 0 , -2 , 8 , -3 , -3 , -1 , -2 , -1 , -2 , -1 , -2 , -2 , 2 , -3
, -1 , -3 , -3 , -3 , -1 , -3 , -3 , -4 , -3 , 4 , 2 , -3 , 1 , 0 , -3 , -2 , -1 , -3 , -1 , 3
, -1 , -2 , -3 , -4 , -1 , -2 , -3 , -4 , -3 , 2 , 4 , -2 , 2 , 0 , -3 , -2 , -1 , -2 , -1 , 1
, -1 , 2 , 0 , -1 , -3 , 1 , 1 , -2 , -1 , -3 , -2 , 5 , -1 , -3 , -1 , 0 , -1 , -3 , -2 , -2
, -1 , -1 , -2 , -3 , -1 , 0 , -2 , -3 , -2 , 1 , 2 , -1 , 5 , 0 , -2 , -1 , -1 , -1 , -1 , 1
, -2 , -3 , -3 , -3 , -2 , -3 , -3 , -3 , -1 , 0 , 0 , -3 , 0 , 6 , -4 , -2 , -2 , 1 , 3 , -1
, -1 , -2 , -2 , -1 , -3 , -1 , -1 , -2 , -2 , -3 , -3 , -1 , -2 , -4 , 7 , -1 , -1 , -4 , -3 , -2
, 1 , -1 , 1 , 0 , -1 , 0 , 0 , 0 , -1 , -2 , -2 , 0 , -1 , -2 , -1 , 4 , 1 , -3 , -2 , -2
, 0 , -1 , 0 , -1 , -1 , -1 , -1 , -2 , -2 , -1 , -1 , -1 , -1 , -2 , -1 , 1 , 5 , -2 , -2 , 0
, -3 , -3 , -4 , -4 , -2 , -2 , -3 , -2 , -2 , -3 , -2 , -3 , -1 , 1 , -4 , -3 , -2 , 11 , 2 , -3
, -2 , -2 , -2 , -3 , -2 , -1 , -2 , -3 , 2 , -1 , -1 , -2 , -1 , 3 , -3 , -2 , -2 , 2 , 7 , -1
, 0 , -3 , -3 , -3 , -1 , -2 , -2 , -3 , -3 , 3 , 1 , -2 , 1 , -1 , -2 , -2 , 0 , -3 , -1 , 4
    ]
    
    HSP={}
    Threshold = 11
    Residue_Length = 3
    
    def __init__(selfparams):
        '''
        Constructor
        '''
        selfparams.Prtn_Alpha_dict={selfparams.Prtn_Alpha[i]:i for i in list(range(len(selfparams.Prtn_Alpha)))}
        selfparams.Prtn_Alpha_length=len(selfparams.Prtn_Alpha)
    
    
    def score(self,resd1,resd2,length=Residue_Length):
        sc=0
        for i in range(length):
            sc+=self.blosum62[
                self.Prtn_Alpha_dict[resd1[i]]+
                self.Prtn_Alpha_dict[resd2[i]]*self.Prtn_Alpha_length
                ]
        return sc
    
    def Try_Gen_Them(self,resd):
        for i in list(itertools.product(self.Prtn_Alpha, repeat=self.Residue_Length-1)):
            yield i[0]+i[1]+resd[2]
            yield resd[0]+i[0]+i[1]
            yield i[0]+resd[1]+i[1]
            
    
    def Generate_Residue_HSPs(self,resd,index):
        if self.HSP.__contains__(resd):
            self.HSP[resd]["Place"].append(index)
            return
        HSP={}
        for i in self.Try_Gen_Them(resd):
            sc=self.score(i, resd)
            if sc > self.Threshold:
                HSP[i]=sc
                HSP["Place"]=[]
                HSP["Place"].append(index)
        self.HSP[resd]=HSP
        
    def Generate_Residue_From_Sequence(self,seq):
        [self.Generate_Residue_HSPs(seq[x:x+self.Residue_Length],x) for x in range(len(seq)) if x+3<=len(seq)]
                