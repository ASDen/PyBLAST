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
    Residue_Length = 3
    Alignment_Count=25


    def __init__(selfparams,HPS,Req):
        '''
        Constructor
        '''
        selfparams.HPS=HPS
        selfparams.Request=Req
        
    def Add_TopHits(self,hit):
        if self.Tfill and hit[3]<self.keys[0]:return
        if hit in self.TopHits:return
        place=bisect(self.keys,hit[3])
        self.TopHits.insert(len(self.keys)- place, hit)
        self.keys.insert(place, hit[3])
        
        if self.Tfill:
            self.TopHits.pop(len(self.keys)-1)
            self.keys.pop(0)
        
        if not self.Tfill and len(self.keys)>=self.Alignment_Count:
            self.Tfill=True
        
    def Extend_Hit(self,seq,st,resd,Sind):
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
            self.Add_TopHits([RsSt,SqSt,SqEnd,Scr,Sind])

    def Get_Record_Hits_with_HPS(self,seq,resd,word,index):
        h=[]
        i=0
        while i < len(seq):
            i = seq.find(resd,i)
            if  i>=0:
                """ Extend hit in both ways """
                self.Extend_Hit(seq, i, word,index)
                i+=self.Residue_Length
            else:
                break
        """ Sort Hits , keep only top scoring ones """
        return h
    
    def GetHits(self,seqs):
        for i in range(len(seqs)):
            if i%100==0:print(i)
            for words in self.HPS:
                for h in self.HPS[words]:
                    if h=="Place": continue
                    self.Get_Record_Hits_with_HPS(seqs[i], h,words,i)
                        
    def Get_Top_Scoring_Alignments(self):
        self.Hits.sort(key=lambda x:x[3],reverse=True)
        i=0
        while i<len(self.Hits)-1:
            if self.Hits[i] == self.Hits[i+1]:
                del self.Hits[i]
            else:
                i=i+1
            
        self.Hits=self.Hits[:self.Alignment_Count]

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
myAnalyzer.GetHits(seqs)
#cProfile.run('myAnalyzer.GetHits(seqs)')
print(str(time.clock())+"  Hits Retrieved")
#myAnalyzer.Get_Top_Scoring_Alignments()
print(myAnalyzer.TopHits)
print(str(time.clock())+" Alignments Ordered "+str(len(myAnalyzer.TopHits)))