#pragma once
class cltm_model
{
public:
	double alpha_phi;
	double alpha_lamda;
	double alpha_gamma;
	double alpha_theta;
	double a;
    int num_topics;
    int num_terms;
	int num_lamda;

	cltm_model(void);
	~cltm_model(void);

    void model_init(double phi,double lamda,double gamma,double theta,int nterm,int ntopics);
};

