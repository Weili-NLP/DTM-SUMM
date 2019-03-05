#coding:utf8

from __future__ import division
import sys
import os
import os.path
import json
import numpy as np
import data
import nltk
from nltk.corpus import PlaintextCorpusReader

def Load_Dict(filename):
    w_dict={}
    with open(filename) as f:
        for line in f:
            data=json.loads(line)
            w_dict[data[0]]=data[1]
    print "Succeed in loading dictionary."        
    return w_dict

def Load_Phil(thetag_file,L,V):
    thetag=np.zeros((L,V))
    l=0
    for topic in file(thetag_file, 'r'):
        thetag[l] = map(float, topic.split())
        l+=1
    return thetag
    
def cal_sentscore(thetag, sent, l):
    sc=0.0
    if len(sent)==0:
        return 0.0
    for w in sent:
        sc+=thetag[l][w]
    return sc/(1.0+0.5*len(sent))

def recover_sent_words(vocab,sent_widx):
    sent=[]
    for w in sent_widx:
        sent.append(vocab[w])
    return sent

def recover_rawtext_sent(raw_text_dict, docid, m):
    s=" "
    content=raw_text_dict[docid]
    return s.join(content[m])

def select_sent(raw_text_dict,vocab, data, thetag, sent_scores, selected_sent_list):
    selnum=20
    for l,lset in enumerate(data.corpus_sent):
        print l,
        for m,doc in enumerate(lset):
            print "  "+str(m)+": ",
            for s,sent in enumerate(doc):
                sent_scores[l][m][s]=cal_sentscore(thetag, sent, l)
            indices = range(len(sent_scores[l][m]))
            indices.sort(lambda x,y: -cmp(sent_scores[l][m][x], sent_scores[l][m][y]))
            for i in range(selnum):
                selected_sent_list[l][m][i]=indices[i]
                print str(indices[i])+" : ",
                #print recover_sent_words(vocab,data.corpus_sent[l][m][indices[i]])
                docid=data.corpus_docids[l][m]
                print recover_rawtext_sent(raw_text_dict, docid, indices[i])
            print '\n'
    return sent_scores, selected_sent_list


def Load_RawText_Collection(root):
    raw_text_dict={}
    corpus_root=root
    wordlists=PlaintextCorpusReader(corpus_root,'.*\.txt')
    porter=nltk.PorterStemmer()
    stopwords=nltk.corpus.stopwords.words('english')
    for each_file in wordlists.fileids():
        content_list=[]
        for each_sent in wordlists.sents(each_file):
            sent_list=[]
            for w in each_sent:
                if len(w)>2 and w.isalpha() and w not in stopwords:
                    sent_list.append(porter.stem(w.lower()))
            if len(sent_list)>0:
                content_list.append(each_sent)
                raw_text_dict[each_file]=content_list
    return raw_text_dict


if __name__ =='__main__':
    docs_fname=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\PaperCollection.txt"
    dict_name=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\Vocab-6\Dictionary.txt"
    cat_root=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\Category"
    #thetag_fname=r"C:\eval\final.thetag"
    dsc_fname=r"C:\eval\final.dsc"
    root=r"C:\Users\hel2\Desktop\F-HLda Model\LDA data set\LDA-27\LDA_Corpus-27"
    selnum=20
    data=data.data()
    data.read_data(dict_name,docs_fname,cat_root)
    w_dict=Load_Dict(dict_name)
    V=len(w_dict)
    vocab=range(V)
    for w, widx in w_dict.items():
        vocab[widx]=w
    #thetag=Load_Phil(thetag_fname,3,V)
    dsc=Load_Phil(dsc_fname,3,V)
    raw_text_dict=Load_RawText_Collection(root)
    sel_sentno_list=[]
    sent_sc=[]
    for l,lset in enumerate(data.corpus_sent):
        l_sc=[]
        l_sel_sent=[]
        for m,doc in enumerate(lset):
            doc_sc=np.zeros(len(doc))
            doc_sel_sent=np.zeros(selnum)
            l_sc.append(doc_sc)
            l_sel_sent.append(doc_sel_sent)
        sel_sentno_list.append(l_sel_sent)
        sent_sc.append(l_sc)

    (sent_sc, sel_sentno_list)=select_sent(raw_text_dict, vocab,data,dsc, sent_sc, sel_sentno_list)
    
    
