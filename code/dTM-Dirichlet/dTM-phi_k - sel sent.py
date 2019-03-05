#coding:utf8
#refer to https://github.com/shuyo/iir/tree/master/lda

import sys
import json
import random
import numpy as np
from time import time
import data

def factorial(base,incm):
    product=1.0
    for i in range(1,incm+1):
        product*=(base+i-1)
    return product

def normalize_p(p):
    sum_p = sum(p)
    for i, pi in enumerate(p):
        p[i] = p[i] / sum_p
    return p

def sample_unnormalize_p(p):
    K=len(p)
    for i in range(K-1):
        p[i+1]+=p[i]
    u=random.random()*p[K-1]
    topic=-1
    for i in range(K):
        if u<p[i]:
            break
    topic=i
    return topic

def SampleFromLogPr(log_pr):
    if len(log_pr)>0:
        log_sum=log_pr[0]
        levels=len(log_pr)
        for i in range(1,levels):
            log_sum=LogSum(log_sum,log_pr[i])
        rand_no=random.random()
        log_exp=math.exp(log_pr[0]-log_sum)
        result=0
        while rand_no>=log_exp:
            result+=1
            log_exp+=math.exp(log_pr[result]-log_sum)
        return result

def LogSum(loga,logb):
    if loga<logb:
        return logb+math.log(1+math.exp(loga-logb))
    else:
        return loga+math.log(1+math.exp(logb-loga))

def gamma(value):
    while value>=1:
        s*=value
        value-=1
    return s

def lngamma(value):
    return math.log(gamma(value-1))
        
class dTM:
    def __init__(self,data,K=20):
        
        
        self._D=data._D
        self._V=data._V
        self._K=K
        self._L=data._L
        self.corpus=data.corpus
        

        self.betab=0.1
        self.betal=0.1
        self.betak=0.1
        self.alpha_lamda=[5.0,3.0,2.0]
        self.alpha_gamma=50/K
        
        self.nb_v=np.zeros(self._V)+self.betab
        self.nl_v=np.zeros((self._L,self._V))+self.betal
        self.nk_v=np.zeros((K,self._V))+self.betak

        self.nb=0+self._V*self.betab
        self.nl=np.zeros(self._L)+self._V*self.betal
        self.nk=np.zeros(K)+self._V*self.betak

        self.nds=[]
        self.ndz=[]
        for l in range(self._L):
            set_ns=[]
            set_nz=[]
            for m in range(len(data.corpus[l])):
                dns=np.zeros(3)+self.alpha_lamda #self.alpha_lamda
                dnz=np.zeros(K)+self.alpha_gamma
                set_ns.append(dns)
                set_nz.append(dnz)
            self.nds.append(set_ns)
            self.ndz.append(set_nz)
            
        self.DWS=[]
        self.DWZ=[]
        for l in range(self._L):
            set_s=[]
            set_z=[]
            for m in range(len(data.corpus[l])):
                doc_s=[]
                doc_z=[]
                for n in range(len(data.corpus[l][m])):
                    widx=data.corpus[l][m][n]
                    s=np.random.randint(0, 3)
                    if s==1:
                        z=np.random.randint(0, K)
                    else:
                        z=-1
                    doc_s.append(s)
                    doc_z.append(z)
                    self.nds[l][m][s]+=1
                    if s==2:
                        self.nb_v[widx]+=1
                        self.nb+=1
                    if s==0:
                        self.nl_v[l][widx]+=1
                        self.nl[l]+=1
                    if s==1:
                        self.nk[z]+=1
                        self.ndz[l][m][z]+=1
                        self.nk_v[z][widx]+=1
                set_s.append(doc_s)
                set_z.append(doc_z)
            self.DWS.append(set_s)
            self.DWZ.append(set_z)
        

    def sample_sz(self,l,m,n): #l,m,n is word index in corpus
        ps=np.zeros(3)
        pz=np.zeros(self._K)
        old_s=self.DWS[l][m][n]
        old_z=self.DWZ[l][m][n]
        widx=self.corpus[l][m][n]
        
        # remove old_s old_z
        self.nds[l][m][old_s]-=1
    
        if old_s==0:
            self.nl[l]-=1
            self.nl_v[l][widx]-=1
        if old_s==1:
            self.nk[old_z]-=1
            self.ndz[l][m][old_z]-=1
            self.nk_v[old_z][widx]-=1
        if old_s==2:
            self.nb-=1
            self.nb_v[widx]-=1

        # cal pz
        for i in range(self._K):
            pz[i]=self.nk_v[i][widx]*self.ndz[l][m][i]/self.nk[i]
        new_z=sample_unnormalize_p(pz)

        #cal ps
        ps[2]=self.nb_v[widx] / self.nb * self.nds[l][m][2] #s=2
        ps[0]=self.nl_v[l][widx] / self.nl[l] * self.nds[l][m][0] #s=0
        ps[1]=self.nk_v[new_z][widx] / self.nk[new_z] * self.nds[l][m][1] #s=1
        new_s=sample_unnormalize_p(ps)
        '''
        print "ps: \n"
        for i in range(len(ps)):
            print str(ps[i])+" "
        print "\n"
        print "pz: \n"
        for i in range(len(pz)):
            print str(pz[i])+" "
        print "\n"
        #dd = raw_input()
        '''
        '''
        psz[0]=self.nb_v[widx]*self.nds[l][m][2]/self.nb #s=2
        psz[1]=self.nl_v[l][widx]*self.nds[l][m][0]/self.nl[l] #s=0
        for i in range(self._K):
            psz[i+2]=self.nk_v[i][widx]*self.nds[l][m][1]*self.ndz[l][m][i]/self.nk[i]

        print "psz: \n"
        for i in range(self._K+2):
            print str(psz[i])+" "
        print "\n"
        dd = raw_input()
        '''

        # add new_sz
        if new_s==2:
            self.DWS[l][m][n]=2
            self.nds[l][m][2]+=1
            self.DWZ[l][m][n]=-1
            self.nb+=1
            self.nb_v[widx]+=1

        elif new_s==0:
            self.DWS[l][m][n]=0
            self.nds[l][m][0]+=1
            self.DWZ[l][m][n]=-1
            self.nl[l]+=1
            self.nl_v[l][widx]+=1
            
        else:
            self.DWS[l][m][n]=1
            self.nds[l][m][1]+=1
            self.DWZ[l][m][n]=new_z
            self.nk[new_z]+=1
            self.ndz[l][m][new_z]+=1
            self.nk_v[new_z][widx]+=1

    def inference(self):
        """learning once iteration over corpus"""
        for l in range(self._L):
            for m in range(len(self.corpus[l])):
                for n in range(len(self.corpus[l][m])):
                    self.sample_sz(l,m,n)

    def topic_worddist(self):
        """get word distribution"""
        return self.nk_v / self.nk[:, np.newaxis]
    def set_worddist(self):
        """get word distribution"""
        return self.nl_v / self.nl[:, np.newaxis]
    def bg_worddist(self):
        """get word distribution"""
        return self.nb_v / self.nb

    def perplexity(self, cps=None):
        if cps == None: cps = self.corpus
        phib = self.bg_worddist()
        phil = self.set_worddist()
        phik = self.topic_worddist()
        log_per = 0
        N = 0
        lamdasum=np.sum(self.alpha_lamda)
        for l, lset in enumerate(cps):
            for m, doc in enumerate(lset):
                d_lamda=self.nds[l][m]/(len(doc)+lamdasum)
                d_gamma=self.ndz[l][m]/(len(doc)+self._K*self.alpha_gamma)
                for n,w in enumerate(doc):
                    pw=0.0
                    pw+=phib[w]*d_lamda[2]
                    pw+=phil[l][w]*d_lamda[0]
                    for k in range(self._K):
                        pw+=phik[k][w]*d_gamma[k]*d_lamda[1]
                    log_per -= np.log(pw)
                N += len(doc)
        return np.exp(log_per / N)

    def save_res(self, dirname, iterN):
        phib_f=dirname+'\\'+str(iterN)+'.phib'
        phil_f=dirname+'\\'+str(iterN)+'.phil'
        phik_f=dirname+'\\'+str(iterN)+'.phik'
        lamda_f=dirname+'\\'+str(iterN)+'.lamda'
        gamma_f=dirname+'\\'+str(iterN)+'.gamma'
        f=open(phib_f,'w')
        for i in range(len(self.nb_v)):
            temp=np.log(self.nb_v[i])-np.log(self.nb)
            if i==0:
                f.write(str(temp))
            else:
                f.write(" "+str(temp))
        f.write("\n")
        f.close()

        f=open(phil_f,'w')
        for i in range(len(self.nl_v)):
            for j in range(len(self.nl_v[i])):
                temp=np.log(self.nl_v[i][j])-np.log(self.nl[i])
                if j==0:
                    f.write(str(temp))
                else:
                    f.write(" "+str(temp))
            f.write("\n")
        f.close()

        f=open(phik_f,'w')
        for i in range(len(self.nk_v)):
            for j in range(len(self.nk_v[i])):
                temp=np.log(self.nk_v[i][j])-np.log(self.nk[i])
                if j==0:
                    f.write(str(temp))
                else:
                    f.write(" "+str(temp))
            f.write("\n")
        f.close()

        f=open(lamda_f,'w')
        for i in range(len(self.nds)):
            for j in range(len(self.nds[i])):
                for s in range(len(self.nds[i][j])):
                    temp=np.log(self.nds[i][j][s])-np.log(len(self.corpus[i][j]))
                    if j==0:
                        f.write(str(temp))
                    else:
                        f.write(" "+str(temp))
                f.write("\n")
        f.close()

        f=open(gamma_f,'w')
        for i in range(len(self.ndz)):
            for j in range(len(self.ndz[i])):
                for s in range(len(self.ndz[i][j])):
                    temp=np.log(self.ndz[i][j][s])-np.log(len(self.corpus[i][j]))
                    if j==0:
                        f.write(str(temp))
                    else:
                        f.write(" "+str(temp))
                f.write("\n")
        f.close()

##############end of class dTM           
    
def dTM_learning(dtm, dirname, iteration):
    pre_perp = dtm.perplexity()
    print ("initial perplexity=%f" % pre_perp)
    for i in range(iteration):
        dtm.inference()
        perp = dtm.perplexity()
        print ("-%d p=%f" % (i + 1, perp))
        if (i%100)==0:
            dtm.save_res(dirname, i)
    dtm.save_res(dirname, iteration)


##############end

                
if __name__ == "__main__":
    docs_fname=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\PaperCollection.txt"
    dict_name=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\Vocab-6\Dictionary.txt"
    cat_root=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\Category"
    data=data.data()
    data.read_data(dict_name,docs_fname,cat_root)
    dirname=r"C:\dTM-gibbs\res"
    dTM = dTM(data,20)
    dTM_learning(dTM, dirname, 3000)
    
            
        

