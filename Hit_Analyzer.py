'''
Created on Jun 22, 2009

@author: Administrator
'''
from bisect import bisect

from QueryProcessor import QueryProcessor
from FASTA_Reader import FASTA_Reader

class Hit_Analyzer(object):
    '''
    Find and Analyze HSP's hits in db Sequences
    '''
    HPS={}
    TopHits=[]
    keys=[]
    Tfill=False
    Request=""
    ReqLen=0
    cSqlLen=0
    Residue_Length = 3
    Alignment_Count=25
    myQP = QueryProcessor()


    def __init__(selfparams,HPS,Req):
        '''
        Constructor
        '''
        selfparams.HPS=HPS
        selfparams.Request=Req
        selfparams.ReqLen=len(Req)-1
        
    def Add_TopHits(self,hit):
        if self.Tfill and hit[3]<self.keys[0]:return
        if hit in self.TopHits:return
        place=bisect(self.keys,hit[3])
        if self.Tfill:
            self.TopHits.insert(self.Alignment_Count - place, hit)
        else:
            self.TopHits.insert(len(self.keys) - place, hit)
        
        self.keys.insert(place, hit[3])
        
        if self.Tfill:
            self.TopHits.pop(self.Alignment_Count-1)
            self.keys.pop(0)
        
        if not self.Tfill and len(self.keys)>=self.Alignment_Count:
            self.Tfill=True
        
    def Extend_Hit(self,seq,st,resd,Sind):
        Ex_Hits=[]
        for RsSt in self.HPS[resd]["Place"]:
            SqSt=st
            SqEnd=st+self.Residue_Length-1
            RsEnd=RsSt+self.Residue_Length-1
            Scr=self.myQP.score(seq[SqSt:SqEnd+1], resd) ##Wasted  recalculations , should be improved
            """ Left Extension """
            while RsSt>0 and SqSt>0:
                nScr=self.myQP.score(seq[SqSt-1], self.Request[RsSt-1], 1)
                if nScr >=0:
                    SqSt-=1
                    RsSt-=1
                    Scr+=nScr
                else:
                    break
            """ Right Extension """
            while RsEnd<self.ReqLen and SqEnd<self.cSqlLen:
                nScr=self.myQP.score(seq[SqEnd+1], self.Request[RsEnd+1], 1)
                if  nScr >=0:
                    SqEnd+=1
                    RsEnd+=1
                    Scr+=nScr
                else:
                    break
            
            self.Add_TopHits([RsSt,SqSt,SqEnd,Scr,Sind])

    def Get_Record_Hits_with_HPS(self,seq,resd,word,index):
        i=0
        while i <= self.cSqlLen:
            i = seq.find(resd,i)
            if  i>=0:
                """ Extend hit in both ways """
                self.Extend_Hit(seq, i, word,index)
                i+=self.Residue_Length
            else:
                break
        
    def GetHits(self,seqs):
        for i in range(len(seqs)):
            if i%100==0:print(i)
            if i>=1000:return
            self.cSqlLen=len(seqs[i])-1
            for words in self.HPS:
                for h in self.HPS[words]:
                    if h=="Place": continue
                    self.Get_Record_Hits_with_HPS(seqs[i], h,words,i)
                        
print("Started")
seqs = []
FR=FASTA_Reader("Mus_musculus.NCBIM30.pep.fa")
#req="AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTACACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTT"
req="HEAAAFLVPVLTHRWNRFAVIVQGEEVTLLMDCEEAAYFMSGLLEEGAGEYDARGYAARTEALAAVVVMDNDSAEVRAYVASADFLDKERAGA"
while True:
    seq=FR.GetNextSequence()
    if seq=='':
        break
    seqs.append(seq)

print("Sequences Loaded")

import time
import cProfile
print(time.clock())
myQP = QueryProcessor()
myQP.Generate_Residue_From_Sequence(req)
#cProfile.run('myQP.Generate_Residue_From_Sequence(req)')
print(str(time.clock())+" Query Processed")
myAnalyzer=Hit_Analyzer(myQP.HSP,req)
#myAnalyzer.GetHits(seqs)
cProfile.run('myAnalyzer.GetHits(seqs)')
print(str(time.clock())+"  Hits Retrieved")
#myAnalyzer.Get_Top_Scoring_Alignments()
print(myAnalyzer.TopHits)
print(str(time.clock())+" Alignments Ordered "+str(len(myAnalyzer.TopHits)))