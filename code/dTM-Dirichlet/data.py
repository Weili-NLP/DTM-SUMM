#coding:utf8

import sys
import os
import os.path
import json

class data:
    def __init__(self):
        self._L=0
        self._V=0
        self._D=0
        self.corpus=[]
        self.corpus_sent=[]
        self.corpus_docids=[] # save the docid for the doc index by l,m

    def load_dict(self,filename):
        w_dict={}
        with open(filename) as f:
            for line in f:
                data=json.loads(line)
                w_dict[data[0]]=data[1]
        print "Succeed in loading dictionary."        
        return w_dict

    def category_dict(self,cat_root):    
        cat_dict={}
        cat_no=0
        for root,dirs,files in os.walk(cat_root):
            for fname in files:
                cat_no+=1
                fname1=fname[:-4]
                fpath=os.path.join(cat_root, fname)
                with open(fpath) as f:
                    for line in f:
                        tag=line[0:3]
                        if tag == '\xef\xbb\xbf':
                            line=line[3:]
                        pos=line.find(':')
                        line=line[0:pos]
                        cat_dict[line]=fname1
        print "Succeed in matching category."        
        return [cat_dict, cat_no]

    def read_data(self,dict_fname,docs_fname,cat_root):
        w_dict=self.load_dict(dict_fname)
        self._V=len(w_dict.keys())
        cat_no=0
        cat_dict={}
        [cat_dict,cat_no]=self.category_dict(cat_root)  #fileid:cat
        self._L=cat_no
        catmap={}
        catmap["LDA_Summarization"]=0
        catmap["LDA_Geo"]=1
        catmap["LDA_Sentiment"]=2
        for i in range(self._L):
            self.corpus.append([])
            self.corpus_sent.append([])
            self.corpus_docids.append([])
        with open(docs_fname) as f:
            dno=0
            for line in f:
                dno+=1
                paper=json.loads(line)
                paper_id=paper['Id']
                l=catmap[cat_dict[paper_id[:-4]]]
                self.corpus_docids[l].append(paper_id)
                paper_content=paper['Content']
                doc=[]
                doc_sent=[]
                for sent in paper_content:
                    sent_word=[]
                    for word in sent:
                        if word in w_dict:
                            w_idx=w_dict[word]
                            doc.append(w_idx)
                            sent_word.append(w_idx)
                    if len(sent_word)!=0:
                        doc_sent.append(sent_word)
                self.corpus[l].append(doc)
                self.corpus_sent[l].append(doc_sent)
            self._D=dno
        print "Succeed in reading corpus."
        print "The number of docs: "+str(self._D)
        print "The number of vobs: "+str(self._V)
        print "The number of cats: "+str(self._L)
        for i in range(self._L):
            print "The number of docs in cat "+str(i)+": "+str(len(self.corpus[i]))

if __name__ =='__main__':
    docs_fname=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\PaperCollection.txt"
    dict_name=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\Vocab-6\Dictionary.txt"
    cat_root=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\Category"
    data=data()
    data.read_data(dict_name,docs_fname,cat_root)
    
