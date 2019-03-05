#pragma once
#include "document.h"
#include "multiset.h"
class corpus
{
public:

	multiset* mlsets;
    int num_terms;  //unique
    int num_docs;   //unique
	int num_multisets;
	int word_total;  // sum of term_counts;
	int doc_total;   // sum of doc_counts;
	int* term_counts;  
	int* doc_counts;

	corpus(void);
	~corpus(void);
};

