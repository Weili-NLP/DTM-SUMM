#pragma once
class suffstats
{
public:
	double** class_word;
    double* class_total;
    double alpha_suffstats;
    int num_docs;

	suffstats(void);
	~suffstats(void);
};

