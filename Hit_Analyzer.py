'''
Created on Jun 22, 2009

@author: Administrator
'''
from QueryProcessor import QueryProcessor
from FASTA_Reader import FASTA_Reader

class Hit_Analyzer(object):
    '''
    Find and Analyze HSP's hits in db Sequences
    '''
    HPS={}
    TopHits=[]
    Request=""
    Residue_Length = 3
    Alignment_Count=25

    def __init__(selfparams,HPS,Req):
        '''
        Constructor
        '''
        selfparams.HPS=HPS
        selfparams.Request=Req
        
    def Extend_Hit(self,seq,st,resd):
        myQP = QueryProcessor()
        Ex_Hits=[]
        for RsSt in self.HPS[resd]["Place"]:
            SqSt=st
            SqEnd=st+self.Residue_Length-1
            RsEnd=RsSt+self.Residue_Length-1
            Scr=myQP.score(seq[SqSt:SqEnd+1], resd) ##Wasted  recalculations , should be improved
            if RsSt>0 and SqSt>0:
                """ Left Extension """
                while RsSt>0 and SqSt>0:
                    nScr=myQP.score(seq[SqSt-1], self.Request[RsSt-1], 1)
                    if nScr >=0:
                        SqSt-=1
                        RsSt-=1
                        Scr+=nScr
                    else:
                        break
            if RsEnd<len(self.Request)-1 and SqEnd<len(seq)-1:
                """ Right Extension """
                while RsEnd<len(self.Request)-1 and SqEnd<len(seq)-1:
                    nScr=myQP.score(seq[SqEnd+1], self.Request[RsEnd+1], 1)
                    if  nScr >=0:
                        SqEnd+=1
                        RsEnd+=1
                        Scr+=nScr
                    else:
                        break
            Ex_Hits.append([RsSt,SqSt,SqEnd,Scr])
        return Ex_Hits
        
    def Get_Record_Hits_with_HPS(self,seq,resd,word,index):
        h=[]
        i=0
        while i < len(seq):
            i = seq.find(resd,i)
            if  i>=0:
                """ Extend hit in both ways """
                for ExHits in self.Extend_Hit(seq, i, word):
                    ExHits.append(index)
                    h.append(ExHits)
                i+=self.Residue_Length
            else:
                break
        """ Sort Hits , keep only top scoring ones """
        return h
    
    def GetHits(self,seqs):
        for i in range(len(seqs)):
            for words in self.HPS:
                for h in self.HPS[words]:
                    if h=="Place": continue
                    hit=self.Get_Record_Hits_with_HPS(seqs[i], h,words,i)
                    if hit != []:
                        self.TopHits.extend(hit)
                        
    def Get_Top_Scoring_Alignments(self):
        self.TopHits.sort(key=lambda x:x[3],reverse=True)
        i=0
        while i<len(self.TopHits)-1:
            if self.TopHits[i] == self.TopHits[i+1]:
                del self.TopHits[i]
            else:
                i=i+1
            
        self.TopHits=self.TopHits[:self.Alignment_Count]

seqs = []
FR=FASTA_Reader("db.fasta")
req="HEAAAFLVPVLTHRWNRFAVIVQGEEVTLLMDCEEAAYFMSGLLEEGAGEYDARGYAARTEALAAVVVMDNDSAEVRAYVASADFLDKERAGA"
while True:
    seq=FR.GetNextSequence()
    if seq=='':
        break
    seqs.append(seq)

import time
import cProfile
print(time.clock())
myQP = QueryProcessor()
#myQP.Generate_Residue_From_Sequence(req)
cProfile.run('myQP.Generate_Residue_From_Sequence(req)')
print(time.clock())
myAnalyzer=Hit_Analyzer(myQP.HSP,req)
print(time.clock())
myAnalyzer.GetHits(seqs)
print(time.clock())
myAnalyzer.Get_Top_Scoring_Alignments()
print(myAnalyzer.TopHits)