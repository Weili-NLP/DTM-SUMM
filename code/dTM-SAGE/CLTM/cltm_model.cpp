#include "StdAfx.h"
#include "cltm_model.h"


cltm_model::cltm_model(void)
{
}


cltm_model::~cltm_model(void)
{
}

void cltm_model::model_init(double phi,double lamda,double gamma,double theta,int nterm,int ntopics)
{
	alpha_phi=phi;
	alpha_lamda=lamda;
	alpha_gamma=gamma;
	alpha_theta=theta;
	a=1;
    num_topics=ntopics;
    num_terms=nterm;
	num_lamda=2;
}