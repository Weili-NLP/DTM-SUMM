#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <time.h>
#include <map>
#include "cltm_model.h"
#include "multiset.h"
#include "corpus.h"
#include "document.h"
#include <gsl/gsl_vector.h>


#define NEWTON_THRESH 1e-3
#define MAX_A_ITER 100
#define MAX_ETA_ITER 100
#define MAX_THETA_ITER 100

#define LAG 5
#define pi  3.1415926535897
#define EM_CONVERGED 1e-15
#define EM_MAX_ITER 10

#define VAR_DOC_CONVERGED 1e-6
#define VAR_SET_CONVERGED 1e-6
#define VAR_MAX_ITER 20
#define VAR_SET_MAX_ITER 4

class estimate
{
public:
	estimate(void);
	~estimate(void);
	bool isnan(double x);
	double sum_array(double* vec, int n);
	double log_sum(double log_a, double log_b);
	/*
	double trigamma(double x);
	double digamma(double x);
	double lgamma(double x);
	double quadgamma(double x);
	*/
	int argmax(double* x, int n);
	double digamma_subtract(double* vec, int index, int n);
	double dl_a(double a, double b, double eta);
	double d2l_a(double a, double b, double eta);
	double opt_a(double init_a, double b, double eta);
	double d_eta(double eta, double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double T,
			   multiset* set,cltm_model* model,int cur_k);
	double d2_eta(double eta, double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			    double T,
			   multiset* set,cltm_model* model,int cur_k);
	double opt_eta(double eta_init,double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			  double T,
			   multiset* set,cltm_model* model,int cur_k);
	//double eta_f(const gsl_vector *eta, void *params);
	//void eta_df(const gsl_vector* eta, void *params, gsl_vector *df);
	//void eta_fdf(const gsl_vector *eta, void* params, double *f, gsl_vector *df);
	double* opt_eta_gsl(double *eta, map<char,void*> *p, int topic_num);
	double* opt_theta_gsl(double *theta, map<char,void*> *p, int topic_num);
	double* opt_a_gsl(double *a, map<char,void*> *p, int topic_num);
	double* opt_b_gsl(double *b, map<char,void*> *p, int topic_num);
	double d_theta(double theta, double**** var_beta_lmn_k,
			   double**** var_lamda_lmn_i,
			   double* T_l,
			   int cur_k,corpus* cps);
	double d2_theta(double theta, double**** var_beta_lmn_k,
			   double**** var_lamda_lmn_i,
			   double* T_l,
			   int cur_k,corpus* cps);
	double opt_theta(double theta, double**** var_beta_lmn_k,
			   double**** var_lamda_lmn_i,
			   double* T_l,
			   int cur_k,corpus* cps);
	double doc_variational_inference(document* doc, cltm_model* model, 
									double* var_alpha_lamda_lm_i, 
									double* var_alpha_gamma_lm_k,
									double** var_beta_lmn_k,
									double** var_lamda_lmn_i,
									double **digamma_subtract_phi_k_v,
									double *dg_topic_k,
									double log_sum_exp_dg_topic);

	double mlset_variational_inference(multiset* set, cltm_model* model, 
									double** var_alpha_lamda_lm_i, 
									double** var_alpha_gamma_lm_k,
									double*** var_beta_lmn_k,
									double*** var_lamda_lmn_i,
									double **digamma_subtract_phi_k_v,
									double *var_alpha_theta_k,
									double *eta_l_k,
									double *var_a_l_k,
									double *var_b_l_k);
	
	double global_variational_inference(corpus* cps, cltm_model* model, 
									double*** var_alpha_lamda_lm_i, 
									double*** var_alpha_gamma_lm_k,
									double**** var_beta_lmn_k,
									double**** var_lamda_lmn_i,
									double **var_alpha_phi_k_v,
									double *var_theta_k,
									double **eta_l_k,
									double **var_a_l_k,
									double **var_b_l_k, char* directory);

   double compute_doc_likelihood(document* doc, cltm_model* model,
					   double* var_alpha_lamda_lm_i, 
					   double* var_alpha_gamma_lm_k,
					   double** var_beta_lmn_k,
					   double** var_lamda_lmn_i,
					   double **digamma_subtract_phi_k_v,
					   double *dg_topic_k,
					  double log_sum_exp_dg_topic);

   double compute_mlset_likelihood(multiset* set, cltm_model* model, 
									double *eta_l_k,
									double *var_a_l_k,
									double *var_b_l_k);

   double compute_global_likelihood(corpus* cps, cltm_model* model, 
									double **var_alpha_phi_k_v,
									double **digamma_subtract_phi_k_v);

   void run_vi(corpus* corpus, char* directory);
};

