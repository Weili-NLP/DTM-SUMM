
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "utils.h"

#include "multiset.h"

#define NEWTON_THRESH 1e-5
#define MAX_ETA_ITER 1000


double d_eta(double eta, double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double a,double b, double T,
			   multiset* set,int cur_k);
double d2_eta(double eta, double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double a,double b, double T,
			   multiset* set,int cur_k);
double opt_eta(double eta_init,double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double a,double b, double T,
			   multiset* set,int cur_k);
//void maximize_alpha(double** gamma, lda_model* model, int num_docs);


