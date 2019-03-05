#include "StdAfx.h"
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <new>
#include "data.h"
data::data(void)
{
}


data::~data(void)
{
}

int* data::copyvec(int* ptr, vector<int> vec){
	ptr=(int*)malloc(vec.size()*sizeof(int));
	for(int i=0;i<vec.size();i++){
		ptr[i]=vec[i];
	}
	return ptr;
}

corpus* data::read_data(string vsm_filename)
{
	vector<int> wordvec, countvec;
	int max_numofsetdocs=50;
    int pos, wid, count,cla_no,item_no;
    corpus* c;
	document* moredocs;
	map<string,int> category_map;
	category_map["LDA_Summarization"]=0;
	category_map["LDA_Geo"]=1;
	category_map["LDA_Sentiment"]=2;

	printf("reading data from %s\n", vsm_filename.c_str());
    c = (corpus*)malloc(sizeof(corpus));
	
	c->num_multisets = 0;
    c->num_terms = 0;
    c->num_docs = 0;
	c->word_total=0;  
	c->doc_total=0; 
   
	// read file
	string line,temp,docid,cla;
	string docnum,vocnum,clanum;
	ifstream ifile;
	ifile.open(vsm_filename.c_str());
	cout<<vsm_filename<<endl;

	if(ifile){
		//read line 1
		getline(ifile,line);
		istringstream instr1(line);
		instr1>>clanum;
		c->num_multisets =stoi(clanum);
		c->mlsets=(multiset*)malloc(c->num_multisets*sizeof(multiset));
		for(int i=0;i<(c->num_multisets);i++){
			c->mlsets[i].size=0;
			//c->mlsets[i].docs=(document*)malloc(max_numofsetdocs * (sizeof(document)));
			c->mlsets[i].docs=new document[max_numofsetdocs];
			if(c->mlsets[i].docs==NULL){
				printf("Error in allocate document pointer!");
			}
		}
		instr1>>docnum;
		// assume doc not overlap in each multiset
		c->num_docs=stoi(docnum);
		c->doc_total=c->num_docs;
		instr1>>vocnum;
		c->num_terms =stoi(vocnum);
		c->term_counts=(int*)malloc(c->num_terms*sizeof(int));  
		for(int i=0;i<c->num_terms;i++){
			c->term_counts[i]=0;
		}
		c->doc_counts=(int*)malloc(c->num_docs*sizeof(int));  
		for(int i=0;i<c->num_docs;i++){
			c->doc_counts[i]=1;
		}
		while(getline(ifile,line)){
				  wordvec.clear();
				  countvec.clear();
                  istringstream instr(line);
				  instr>>docid;
				  instr>>cla;
				  cla_no=category_map[cla];
				  item_no=c->mlsets[cla_no].size;
				  if(item_no==max_numofsetdocs){
					  max_numofsetdocs=2*max_numofsetdocs;
					  moredocs=(document*)realloc(c->mlsets[cla_no].docs,max_numofsetdocs*sizeof(document));
					  if(moredocs!=NULL){
						c->mlsets[cla_no].docs=moredocs;
					  }
				  }
				  
				  
				  c->mlsets[cla_no].docs[item_no].length=0;
				  c->mlsets[cla_no].docs[item_no].total=0;
				  c->mlsets[cla_no].docs[item_no].docid=docid;
				  c->mlsets[cla_no].docs[item_no].cla=cla;
				  
				  while(instr>>temp){
					  pos=temp.find(':');
					  wid=stoi(temp.substr(0,pos));
					  count=stoi(temp.substr(pos+1));
					  c->word_total+=count;
					  c->term_counts[wid]+=count;
					  c->mlsets[cla_no].docs[item_no].total+=count;
					  c->mlsets[cla_no].docs[item_no].length++;
					  wordvec.push_back(wid);
					  countvec.push_back(count);
				  }
				  c->mlsets[cla_no].docs[item_no].words = (int*)malloc(wordvec.size()*sizeof(int));
				  c->mlsets[cla_no].docs[item_no].counts = (int*)malloc(countvec.size()*sizeof(int));
				  for(int i=0;i<wordvec.size();i++){
						c->mlsets[cla_no].docs[item_no].words[i] = wordvec[i];
						c->mlsets[cla_no].docs[item_no].counts[i] = countvec[i];
				  }
				  c->mlsets[cla_no].size++;
		}
	}else{
		cout<<"File Not Found!"<<endl;
	}
	ifile.close();
	printf("number of docs: %d\n", c->num_docs);
	printf("number of terms: %d\n", c->num_terms);

	/*
	for (map<string,int>::iterator it=category_map.begin(); it!=category_map.end(); ++it){
		printf("number of docs in multiset %d: %d\n", it->second, c->mlsets[it->second].size);
		for(int i=0;i<c->mlsets[it->second].size;i++){
			printf("number of unque words in doc %d multiset %d: %d\n", i, it->second, c->mlsets[it->second].docs[i].length);
			if(i==0){
				for(int j=0;j<c->mlsets[it->second].docs[i].length;j++){
					printf("number of word %d : %d\n", c->mlsets[it->second].docs[i].words[j],c->mlsets[it->second].docs[i].counts[j]);
				}
			}
		}
	}
	*/
    return(c);
}
