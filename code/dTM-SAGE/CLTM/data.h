#pragma once
#include "corpus.h"
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
class data
{
public:
	
	data(void);
	~data(void);
	int* data::copyvec(int* ptr, vector<int> vec);
	corpus* read_data(string data_filename);
};

