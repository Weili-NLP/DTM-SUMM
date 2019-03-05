#pragma once
#include "word.h"
#include <vector>

using namespace std;

class document
{
public:
	string docid;
	string cla;
	int* words;
	int* counts;
	int length;
	int total;
 
    
	document(void);
	~document(void);
};

