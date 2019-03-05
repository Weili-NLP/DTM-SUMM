#include "StdAfx.h"
#include <iostream>
#include <random>
#include "estimate.h"
//#include "utils.h"
//#include "opt_a.h"
//#include "opt_eta.h"
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_multimin.h>


using namespace std;

estimate::estimate(void)
{
}


estimate::~estimate(void)
{
}



 /**
   * Proc to calculate the value of the trigamma, the second
   * derivative of the loggamma function. Accepts positive matrices.
   * From Abromowitz and Stegun.  Uses formulas 6.4.11 and 6.4.12 with
   * recurrence formula 6.4.6.  Each requires workspace at least 5
   * times the size of X.
   *
   **/

double trigamma(double x)
{
    double p;
    int i;

    x=x+6;
    p=1/(x*x);
    p=(((((0.075757575757576*p-0.033333333333333)*p+0.0238095238095238)
         *p-0.033333333333333)*p+0.166666666666667)*p+1)/x+0.5*p;
    for (i=0; i<6 ;i++)
    {
        x=x-1;
        p=1/(x*x)+p;
    }
    return(p);
}

double quadgamma(double x)
{
	return gsl_sf_psi_n(2,x);
	//return 0;
}
/*
 * taylor approximation of first derivative of the log gamma function
 *
 */

double digamma(double x)
{
    double p;
    x=x+6;
    p=1/(x*x);
    p=(((0.004166666666667*p-0.003968253986254)*p+
	0.008333333333333)*p-0.083333333333333)*p;
    p=p+log(x)-0.5/x-1/(x-1)-1/(x-2)-1/(x-3)-1/(x-4)-1/(x-5)-1/(x-6);
    return p;
}


double lgamma(double x)
{
     double z=1/(x*x);

    x=x+6;
    z=(((-0.000595238095238*z+0.000793650793651)
	*z-0.002777777777778)*z+0.083333333333333)/x;
    z=(x-0.5)*log(x)-x+0.918938533204673+z-log(x-1)-
	log(x-2)-log(x-3)-log(x-4)-log(x-5)-log(x-6);
    return z;
}



///// model output

void save_theta0(char* filename, double* var_theta_k, int num_topics)
{
    FILE* fileptr;
    int k;
    fileptr = fopen(filename, "w");

   
	fprintf(fileptr, "%5.10f",var_theta_k[0]);
	for (k = 1; k < num_topics; k++)
	{
	    fprintf(fileptr, " %5.10f", var_theta_k[k]);
	}
	fprintf(fileptr, "\n");
    
    fclose(fileptr);
}

void save_thetag(char* filename, double* var_theta_k, double** eta_l_k,int num_set, int num_topics)
{
    FILE* fileptr;
    int d,k;
    fileptr = fopen(filename, "w");

    for (d = 0; d < num_set; d++)
    {
	fprintf(fileptr, "%5.10f", exp(eta_l_k[d][0]+var_theta_k[0]));
	for (k = 1; k < num_topics; k++)
	{
	    fprintf(fileptr, " %5.10f", exp(eta_l_k[d][k]+var_theta_k[k]));
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);
}


void save_alpha_phi(char* filename, double**var_alpha_phi_k_v, int num_topics, int num_terms)
{
    FILE* fileptr;
    int k,v;
    fileptr = fopen(filename, "w");

    for (k = 0; k < num_topics; k++)
    {
	fprintf(fileptr, "%5.10f", var_alpha_phi_k_v[k][0]);
	for (v = 1; v < num_terms; v++)
	{
	    fprintf(fileptr, " %5.10f", var_alpha_phi_k_v[k][v]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);
}

//////end of model output

//gsl_opt_b
double b_f(const gsl_vector *b, void *params){
	map<char,void*> *p=(map<char,void*>*)params;
	double* eta = (double*) p->at('a');
	double* a= (double*) p->at('b');
	int num_topics=b->size;
	int k;
	double sum=0.0;
	double r=1.0;
	double x=0.0;
	for(k=0;k<num_topics;k++){
		x=gsl_vector_get(b,k);
		sum+=0.5*log(x)-0.5*pow(eta[k],2)*(1/(x*(a[k]-1)))-r*x*a[k];
	}
	return -sum;
}

void b_df(const gsl_vector* b, void *params, gsl_vector *df){
	map<char,void*> *p=(map<char,void*>*)params;
	double* eta = (double*) p->at('a');
	double* a= (double*) p->at('b');
	int num_topics=b->size;
	int k;
	double sum=0.0;
	double r=1.0;
	double x=0.0;
	for(k=0;k<num_topics;k++){
		x=gsl_vector_get(b,k);
		sum = 0.5/x-r*a[k]+0.5*pow(eta[k],2)*(1/(a[k]-1))*(1/pow(x,2));
		gsl_vector_set(df,k,-sum);
	}
	
}
void b_fdf(const gsl_vector *b, void* params, double *f, gsl_vector *df){
	*f=b_f(b,params);
	b_df(b,params,df);
}

double* estimate::opt_b_gsl(double *b, map<char,void*> *p, int topic_num){
	int status;
	int iter=0;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	
	gsl_vector *x;
	gsl_multimin_function_fdf b_func;

	b_func.n=topic_num;
	b_func.f=&b_f;
	b_func.df=&b_df;
	b_func.fdf=&b_fdf;
	b_func.params=(void*)p;

	x=gsl_vector_alloc(topic_num);
	for(int k=0;k<topic_num;k++){
		gsl_vector_set(x,k,b[k]);
	}

	//T=gsl_multimin_fdfminimizer_conjugate_fr;
	//T=gsl_multimin_fdfminimizer_conjugate_pr;
	T=gsl_multimin_fdfminimizer_vector_bfgs;

	s=gsl_multimin_fdfminimizer_alloc(T,topic_num);
	gsl_multimin_fdfminimizer_set(s,&b_func,x,0.01,1e-3);
	
	do{
		iter++;
		status=gsl_multimin_fdfminimizer_iterate(s);
		if(status)
			break;
		status=gsl_multimin_test_gradient(s->gradient,NEWTON_THRESH);
	}while((status==GSL_CONTINUE) && (iter < MAX_A_ITER));

	for(int k=0;k<topic_num;k++){
		b[k]=gsl_vector_get(s->x,k);
	}
	return b;
}

//end of gsl_opt_b





///////

//gsl_opt_a
double a_f(const gsl_vector *a, void *params){
	map<char,void*> *p=(map<char,void*>*)params;
	double* eta = (double*) p->at('a');
	double* b= (double*) p->at('b');
	int num_topics=a->size;
	int k;
	double sum=0.0;
	double r=1.0;
	double x=0.0;
	for(k=0;k<num_topics;k++){
		x=gsl_vector_get(a,k);
		sum+=(0.5-x)*digamma(x)-0.5*pow(eta[k],2)*(1/(b[k]*(x-1)))-r*x*b[k]+x+lgamma(x);
	}
	return -sum;
}

void a_df(const gsl_vector* a, void *params, gsl_vector *df){
	map<char,void*> *p=(map<char,void*>*)params;
	double* eta = (double*) p->at('a');
	double* b= (double*) p->at('b');
	int num_topics=a->size;
	int k;
	double sum=0.0;
	double r=1.0;
	double x=0.0;
	for(k=0;k<num_topics;k++){
		x=gsl_vector_get(a,k);
		sum=(0.5-x)*trigamma(x)+0.5*pow(eta[k],2)*(1/b[k])*(1/pow((x-1),2))-r*b[k]+1;
		gsl_vector_set(df,k,-sum);
	}
	
}
void a_fdf(const gsl_vector *a, void* params, double *f, gsl_vector *df){
	*f=a_f(a,params);
	a_df(a,params,df);
}

double* estimate::opt_a_gsl(double *a, map<char,void*> *p, int topic_num){
	int status;
	int iter=0;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	
	gsl_vector *x;
	gsl_multimin_function_fdf a_func;

	a_func.n=topic_num;
	a_func.f=&a_f;
	a_func.df=&a_df;
	a_func.fdf=&a_fdf;
	a_func.params=(void*)p;

	x=gsl_vector_alloc(topic_num);
	for(int k=0;k<topic_num;k++){
		gsl_vector_set(x,k,a[k]);
	}

	//T=gsl_multimin_fdfminimizer_conjugate_fr;
	//T=gsl_multimin_fdfminimizer_conjugate_pr;
	T=gsl_multimin_fdfminimizer_vector_bfgs;

	s=gsl_multimin_fdfminimizer_alloc(T,topic_num);
	gsl_multimin_fdfminimizer_set(s,&a_func,x,0.01,1e-3);
	
	do{
		iter++;
		status=gsl_multimin_fdfminimizer_iterate(s);
		if(status)
			break;
		status=gsl_multimin_test_gradient(s->gradient,NEWTON_THRESH);
	}while((status==GSL_CONTINUE) && (iter < MAX_A_ITER));

	for(int k=0;k<topic_num;k++){
		a[k]=gsl_vector_get(s->x,k);
	}
	return a;
}

//end of gsl_opt_a

//////
//gsl_opt_theta
double theta_f(const gsl_vector *theta, void *params){
	map<char,void*> *p=(map<char,void*>*)params;
	double** eta = (double**) p->at('a');
	double**** lamda= (double****) p->at('b');
	double**** beta = (double****) p->at('c');
	corpus* cps = (corpus*) p->at('d');
	double** T;
	double* log_sum_exp_l;
	int l,m,n,k;
	double res=0.0, sum=0.0, log_sum_exp=0.0;
	int num_topics=theta->size;
	T=(double**)malloc(cps->num_multisets * sizeof(double*));
	log_sum_exp_l=(double*)malloc(cps->num_multisets * sizeof(double));
	for(l=0;l<cps->num_multisets;l++){
		T[l]=(double*)malloc(num_topics * sizeof(double));
	}
	for(l=0;l<cps->num_multisets;l++){
		log_sum_exp_l[l]=0.0;
		for(k=0;k<num_topics;k++){
			log_sum_exp_l[l]+= exp(gsl_vector_get(theta,k)+eta[l][k]);
		}
		log_sum_exp_l[l]=log(log_sum_exp_l[l]);
	}
	for(l=0;l<cps->num_multisets;l++){
		for(k=0;k<num_topics;k++){
			T[l][k]=gsl_vector_get(theta,k)-log_sum_exp_l[l];
		}
	}
	res=0.0;
	for(l=0;l<cps->num_multisets;l++){
		for(m=0;m<cps->mlsets[l].size;m++){
			for(n=0;n<cps->mlsets[l].docs[m].length;n++){
				sum=0.0;
				for(k=0;k<num_topics;k++){
					sum+=beta[l][m][n][k]*T[l][k];
				}
				res+=lamda[l][m][n][0]*sum;
			}
		}
	}

	for(l=0;l<cps->num_multisets;l++){
		free (T[l]);
	}
	free (T);
	free (log_sum_exp_l);
	return -res;
}

void theta_df(const gsl_vector* theta, void *params, gsl_vector *df){
	map<char,void*> *p=(map<char,void*>*)params;
	double** eta = (double**) p->at('a');
	double**** lamda= (double****) p->at('b');
	double**** beta = (double****) p->at('c');
	corpus* cps = (corpus*) p->at('d');
	double** T;
	double* sum_exp_l;
	
	int l,m,n,k,ck;
	double res=0.0, sum=0.0, sum_exp=0.0;
	int num_topics=theta->size;

	T=(double**)malloc(cps->num_multisets * sizeof(double*));
	sum_exp_l=(double*)malloc(cps->num_multisets * sizeof(double));
	for(l=0;l<cps->num_multisets;l++){
		T[l]=(double*)malloc(num_topics * sizeof(double));
	}
	for(l=0;l<cps->num_multisets;l++){
		sum_exp_l[l]=0.0;
		for(k=0;k<num_topics;k++){
			sum_exp_l[l]+= exp(gsl_vector_get(theta,k)+eta[l][k]);
		}
	}
	for(l=0;l<cps->num_multisets;l++){
		for(k=0;k<num_topics;k++){
			T[l][k]=exp(gsl_vector_get(theta,k)+eta[l][k])/sum_exp_l[l];
		}
	}

	for(ck=0;ck<num_topics;ck++){
		res=0.0;
		for(l=0;l<cps->num_multisets;l++){
			for(m=0;m<cps->mlsets[l].size;m++){
				for(n=0;n<cps->mlsets[l].docs[m].length;n++){
					sum=0.0;
					for(k=0;k<num_topics;k++){
						if(ck==k){
							sum+=beta[l][m][n][k]*(1-T[l][ck]);
						}else{
							sum-=beta[l][m][n][k]*T[l][ck];
						}
					}
					res+=lamda[l][m][n][0]*sum;
				}
			}
		}
		gsl_vector_set(df,ck,-res);
	}

	for(l=0;l<cps->num_multisets;l++){
		free (T[l]);
	}
	free (T);
	free (sum_exp_l);
}
void theta_fdf(const gsl_vector *theta, void* params, double *f, gsl_vector *df){
	*f=theta_f(theta,params);
	theta_df(theta,params,df);
}

double* estimate::opt_theta_gsl(double *theta, map<char,void*> *p, int topic_num){
	int status;
	int iter=0;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	
	gsl_vector *x;
	gsl_multimin_function_fdf theta_func;

	theta_func.n=topic_num;
	theta_func.f=&theta_f;
	theta_func.df=&theta_df;
	theta_func.fdf=&theta_fdf;
	theta_func.params=(void*)p;

	x=gsl_vector_alloc(topic_num);
	for(int k=0;k<topic_num;k++){
		gsl_vector_set(x,k,theta[k]);
	}

	//T=gsl_multimin_fdfminimizer_conjugate_fr;
	//T=gsl_multimin_fdfminimizer_conjugate_pr;
	T=gsl_multimin_fdfminimizer_vector_bfgs;
	s=gsl_multimin_fdfminimizer_alloc(T,topic_num);
	gsl_multimin_fdfminimizer_set(s,&theta_func,x,0.01,1e-3);
	
	do{
		iter++;
		status=gsl_multimin_fdfminimizer_iterate(s);
		if(status)
			break;
		status=gsl_multimin_test_gradient(s->gradient,NEWTON_THRESH);
	}while((status==GSL_CONTINUE) && (iter < MAX_THETA_ITER));

	for(int k=0;k<topic_num;k++){
		theta[k]=gsl_vector_get(s->x,k);
	}
	return theta;
}

//end of gsl_opt_theta

///////

//gsl_opt_eta
double eta_f(const gsl_vector *eta, void *params){
	map<char,void*> *p=(map<char,void*>*)params;
	double* theta = (double*) p->at('a');
	double*** lamda= (double***) p->at('b');
	multiset* set = (multiset*) p->at('d');
	double* T;
	int m,n,v;
	double res=0.0, sum=0.0, log_sum_exp=0.0;
	int num=eta->size;
	T=(double*)malloc(num*sizeof(double));
	for(v=0;v<num;v++){
		log_sum_exp+=exp(theta[v]+gsl_vector_get(eta,v));
	}
	log_sum_exp=log(log_sum_exp);
	for(v=0;v<num;v++){
		T[v]=gsl_vector_get(eta,v)-log_sum_exp;
	}
	res=0.0;
	for(m=0;m<set->size;m++){
		for(n=0;n<set->docs[m].length;n++){
			res += (set->docs[m].counts[n])*lamda[m][n][0]*T[set->docs[m].words[n]];
		}
	}
	for(v=0;v<num;v++){
		res-=0.5*pow(gsl_vector_get(eta,v),4);
	}
	free (T);
	return -res;
}

void eta_df(const gsl_vector* eta, void *params, gsl_vector *df){
	map<char,void*> *p=(map<char,void*>*)params;
	double* theta = (double*) p->at('a');
	double*** lamda= (double***) p->at('b');
	//double*** beta = (double***) p->at('c');
	multiset* set = (multiset*) p->at('d');
	double* T;
	int m,n,v,cv;
	double res=0.0, sum=0.0, sum_exp=0.0;
	int num=eta->size;
	T=(double*)malloc(num*sizeof(double));
	for(v=0;v<num;v++){
		sum_exp+=exp(theta[v]+gsl_vector_get(eta,v));
	}
	for(v=0;v<num;v++){
		T[v]=exp(theta[v]+gsl_vector_get(eta,v))/sum_exp;
	}
	for(cv=0;cv<num;cv++){
		res=0.0;
		for(m=0;m<set->size;m++){
			for(n=0;n<set->docs[m].length;n++){
				if(set->docs[m].words[n]==cv){
					res += (set->docs[m].counts[n])*lamda[m][n][0]*(1-T[cv]);
				}else{
					res -= (set->docs[m].counts[n])*lamda[m][n][0]*T[cv];
				}
			}
		}
		res-=2*pow(gsl_vector_get(eta,cv),3);
		gsl_vector_set(df,cv,-res);
	}

	free (T);
}
void eta_fdf(const gsl_vector *eta, void* params, double *f, gsl_vector *df){
	*f=eta_f(eta,params);
	eta_df(eta,params,df);
}

double* estimate::opt_eta_gsl(double *eta, map<char,void*> *p, int num){
	int status;
	int iter=0;
	const gsl_multimin_fdfminimizer_type *T;
	gsl_multimin_fdfminimizer *s;
	
	gsl_vector *x;
	gsl_multimin_function_fdf eta_func;

	eta_func.n=num;
	eta_func.f=&eta_f;
	eta_func.df=&eta_df;
	eta_func.fdf=&eta_fdf;
	eta_func.params=(void*)p;

	x=gsl_vector_alloc(num);
	for(int k=0;k<num;k++){
		gsl_vector_set(x,k,eta[k]);
	}

	//T=gsl_multimin_fdfminimizer_conjugate_fr;
	//T=gsl_multimin_fdfminimizer_conjugate_pr;
	T=gsl_multimin_fdfminimizer_vector_bfgs;
	s=gsl_multimin_fdfminimizer_alloc(T,num);
	gsl_multimin_fdfminimizer_set(s,&eta_func,x,0.01,1e-3);
	
	do{
		iter++;
		status=gsl_multimin_fdfminimizer_iterate(s);
		if(status)
			break;
		status=gsl_multimin_test_gradient(s->gradient,NEWTON_THRESH);
	}while((status==GSL_CONTINUE) && (iter < MAX_ETA_ITER));

	for(int k=0;k<num;k++){
		eta[k]=gsl_vector_get(s->x,k);
	}
	return eta;
}

//end of gsl_opt_eta

//util

bool estimate::isnan(double x) { return x != x; }

double estimate::sum_array(double* vec, int n){
	double sum=0.0;
	for(int i=0;i<n;i++)
		sum=sum+vec[i];
	return (sum);
}
/*
 * given log(a) and log(b), return log(a + b)
 *
 */

double estimate::log_sum(double log_a, double log_b)
{
  double v;

  if (log_a < log_b)
  {
      v = log_b+log(1 + exp(log_a-log_b));
  }
  else
  {
      v = log_a+log(1 + exp(log_b-log_a));
  }
  return(v);
}

// give a_1, ..., a_n,


/*
 * argmax
 *
 */

int estimate::argmax(double* x, int n)
{
    int i;
    double max = x[0];
    int argmax = 0;
    for (i = 1; i < n; i++)
    {
        if (x[i] > max)
        {
            max = x[i];
            argmax = i;
        }
    }
    return(argmax);
}


//end of util

// opt_a

double estimate::dl_a(double a, double b, double eta)
{ 
	double r=1.0;
	return (0.5+a)*trigamma(a) - 0.5*pow(eta,2)*(1/b)*(1/pow((a-1),2)) - 1;
	//return (0.5-a)*trigamma(a) + 0.5*pow(eta,2)*(1/b)*(1/pow((a-1),2)) - r*b + 1;
	//return (0.5-a)*trigamma(a)- r*b + 1;
}

double estimate::d2l_a(double a, double b, double eta)
{ 
	return (0.5+a)*quadgamma(a)+ pow(eta,2)*(1/b)*(1/pow((a-1),3));
	//return (0.5-a)*quadgamma(a)- trigamma(a) - pow(eta,2)*(1/b)*(1/pow((a-1),3));
	//return (0.5-a)*quadgamma(a)- trigamma(a) - b;
}


/*
 * newtons method 
 * ensure a not equals 1.0
 */

double estimate::opt_a(double init_a, double b, double eta)
{
    double a;
    double df, d2f;
    int iter = 0;
	init_a=10;

    a = init_a;
    do
    {
        iter++;
        //a = exp(log_a);
        if (isnan(a)|| (a<=1.0))
        {
            init_a = init_a * 2;
            //printf("warning : a is nan; new init = %5.5f\n", init_a);
            a = init_a;
            //log_a = log(a);
        }
        df = dl_a(a, b, eta);
        d2f = d2l_a(a, b, eta);
        a = a - df/d2f;
        //printf("alpha maximization : %5.5f   %5.5f\n", f, df);
		
    }
    while ( (((fabs(df) < NEWTON_THRESH) && (a>1.0))==false) && (((iter > MAX_A_ITER) && (a>1.0))==false));
    return(a);
}
//end of opt_a

//opt_theta
/*
var_theta_k[k],  T_l, var_beta_lmn_k,
			   var_lamda_lmn_i,k,cps
*/
double estimate::d_theta(double theta, double**** var_beta_lmn_k,
			   double**** var_lamda_lmn_i,
			   double* T_l,
			   int cur_k,corpus* cps)
{ 
	int l,m,n;

	double sum,res;
	res=0.0;
	for(l = 0; l < cps->num_multisets; l++){
		sum=0.0;
		for(m = 0; m < cps->mlsets[l].size; m++){
			for(n=0;n< cps->mlsets[l].docs[m].length;n++){
				sum+=(var_lamda_lmn_i[l][m][n][0]*var_beta_lmn_k[l][m][n][cur_k]);
			}
		}
		res+=sum*(1-T_l[l]);
	}
	
	return res;
}

double estimate::d2_theta(double theta, double**** var_beta_lmn_k,
			   double**** var_lamda_lmn_i,
			   double* T_l,
			   int cur_k,corpus* cps)
{  int l,m,n;

	double sum,res;
	res=0.0;
	for(l = 0; l < cps->num_multisets; l++){
		sum=0.0;
		for(m = 0; m < cps->mlsets[l].size; m++){
			for(n=0;n< cps->mlsets[l].docs[m].length;n++){
				sum+=(var_lamda_lmn_i[l][m][n][0]*var_beta_lmn_k[l][m][n][cur_k]);
			}
		}
		res+=sum*T_l[l]*(T_l[l]-1);
	}
	
	return res;
}

/*
 * newtons method
 *
 */
//sth wrong with this theta derivative
double estimate::opt_theta(double theta_init, double**** var_beta_lmn_k,
			   double**** var_lamda_lmn_i,
			   double* T_l,
			   int cur_k,corpus* cps)
{

	double theta, log_theta;
    double  df, d2f;
    int iter = 0;

    log_theta = log(theta_init);
    do
    {
        iter++;
        theta = exp(log_theta);
        if (isnan(theta))
        {
            theta_init = theta_init * 2;
            printf("warning : theta is nan; new init = %5.5f\n", theta_init);
            theta = theta_init;
            log_theta = log(theta);
        }
        df = d_theta(theta, var_beta_lmn_k,
			   var_lamda_lmn_i,T_l,cur_k,cps);
        d2f = d2_theta(theta, var_beta_lmn_k,
			   var_lamda_lmn_i,T_l,cur_k,cps);
        log_theta = log_theta - df/(d2f * theta + df);
       // printf("theta maximization :  %5.5f\n", df);
    }
    while ((fabs(df) > NEWTON_THRESH) && (iter < MAX_THETA_ITER));
    return(exp(log_theta));
}

//end of opt_theta


//opt_eta
double estimate::d_eta(double eta, double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double T,
			   multiset* set,cltm_model* model,int cur_k)
{ 
	int k,m,n;
	double sum=0.0;
	double res=0.0;
	for(m = 0; m < set->size; m++){
		for(n=0;n<set->docs[m].length;n++){
			sum=0.0;
			for(k = 0; k < model->num_topics; k++){
				if(k==cur_k){
					sum+=var_beta_lmn_k[m][n][k]*(1-T);
				}else{
					sum-=var_beta_lmn_k[m][n][k]*T;
				}
			}			
			res+=sum*var_lamda_lmn_i[m][n][0];
		}
	}
	return res-2*pow(eta,3);
}

//sth wrong with this 2nd derivative of eta
double estimate::d2_eta(double eta, double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double T,
			   multiset* set,cltm_model* model,int cur_k)
{   int m,n;
	double res=0.0;
	for(m = 0; m < set->size; m++){
		for(n=0;n<set->docs[m].length;n++){
			res+=var_lamda_lmn_i[m][n][0]*var_beta_lmn_k[m][n][cur_k];
		}
	}
	return res*(T*(T-1))-6*pow(eta,2);
}

/*
 * newtons method
 *
 */

double estimate::opt_eta(double eta_init,double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double T,
			   multiset* set,cltm_model* model,int cur_k)
{
	
    double eta;
    double df, d2f;
    int iter = 0;
	
	eta_init=10;

    //log_eta = log(eta_init);
	eta=eta_init;
    do
    {
        iter++;
        //eta = exp(log_eta);
        if (isnan(eta)|| (eta==0.0))
        {
            eta_init = (eta_init+0.1)*2;
            printf("warning : eta is nan; new init = %5.5f\n", eta_init);
            eta =eta_init;
            //log_eta = log(eta);
        }
   
        df = d_eta(eta, var_beta_lmn_k,
			   var_lamda_lmn_i,
			  T,set,model,cur_k);
        d2f = d2_eta( eta,  var_beta_lmn_k,
			    var_lamda_lmn_i,
			  T,set,model,cur_k);

        eta = eta - df/d2f ;
       // printf("alpha maximization : %5.5f   %5.5f\n", f, df);
    }
    while ((((fabs(df) < NEWTON_THRESH) && (eta!=0.0))==false) && (((iter > MAX_ETA_ITER) && (eta!=0.0))==false));
    return(eta);
}

//end of eta



double estimate::digamma_subtract(double* vec, int index, int n){
	return digamma(vec[index])-digamma(sum_array(vec,n));
}


/*
 * compute doc-ELOB
 *
 */

double estimate::compute_doc_likelihood(document* doc, cltm_model* model,
					   double* var_alpha_lamda_lm_i, 
					   double* var_alpha_gamma_lm_k,
					   double** var_beta_lmn_k,
					   double** var_lamda_lmn_i,
					   double **digamma_subtract_phi_k_v,
					   double *dg_vob_v,
					  double log_sum_exp_dg_vob)
{
    double likelihood = 0, dig_lamda_sum = 0, dig_gammma_sum = 0, var_gamma_sum = 0, var_lamda_sum = 0;
	double sum=0.0, sum1=0.0;
	double* dig_lamda=(double*)malloc(sizeof(double)*(model->num_lamda));
	double* dig_gamma=(double*)malloc(sizeof(double)*(model->num_topics));
    int i, k, n;

    for (k = 0; k < model->num_topics; k++)
    {
		dig_gamma[k] = digamma(var_alpha_gamma_lm_k[k]);
		var_gamma_sum += var_alpha_gamma_lm_k[k];
    }
    dig_gammma_sum = digamma(var_gamma_sum);

	likelihood =
		lgamma(model->alpha_gamma * model -> num_topics)
		- model -> num_topics * lgamma(model->alpha_gamma)
		- (lgamma(var_gamma_sum));

	for (i = 0; i < model->num_lamda; i++)
    {
		dig_lamda[i] = digamma(var_alpha_lamda_lm_i[i]);
		var_lamda_sum += var_alpha_lamda_lm_i[i];
    }
    dig_lamda_sum = digamma(var_lamda_sum);

    likelihood +=
		(lgamma(model->alpha_lamda * model -> num_lamda)
		- model -> num_lamda * lgamma(model->alpha_lamda)
		- (lgamma(var_lamda_sum)));

    for (k = 0; k < model->num_topics; k++)
    {
		likelihood +=
			((model->alpha_gamma - 1)*(dig_gamma[k] - dig_gammma_sum) + lgamma(var_alpha_gamma_lm_k[k])
			 - (var_alpha_gamma_lm_k[k] - 1)*(dig_gamma[k] - dig_gammma_sum));
	}
	for (i = 0; i < model->num_lamda; i++)
    {
		likelihood +=
			((model->alpha_lamda - 1)*(dig_lamda[i] - dig_lamda_sum) + lgamma(var_alpha_lamda_lm_i[i])
			 - (var_alpha_lamda_lm_i[i] - 1)*(dig_lamda[i] - dig_lamda_sum));
	}


	for (n = 0; n < doc->length; n++)
	{
		for(i = 0; i < model->num_lamda; i++){
			 likelihood += (doc->counts[n]*
                           var_lamda_lmn_i[n][i] * (dig_lamda[i] - dig_lamda_sum - log(var_lamda_lmn_i[n][i])));
		}
		for (k = 0; k < model->num_topics; k++)
		{
			 likelihood += (doc->counts[n]*
						   var_beta_lmn_k[n][k] * (dig_gamma[k] - dig_gammma_sum- log(var_beta_lmn_k[n][k])));
		}
		for(i = 0; i < model->num_lamda; i++){
			if(i==1){
				sum=0.0;
				for(k = 0; k < model->num_topics; k++){
					sum += (var_beta_lmn_k[n][k]* (digamma_subtract_phi_k_v[k][doc->words[n]]));
				}
				likelihood += (doc->counts[n]* var_lamda_lmn_i[n][1]* sum);
			}
			else if(i==0){
				likelihood += (doc->counts[n]* var_lamda_lmn_i[n][0]* (dg_vob_v[doc->words[n]]-log_sum_exp_dg_vob));
			}
		}
    }
   
	free (dig_lamda);
	free (dig_gamma);
    return(likelihood);
}


/*
 * compute mlset-ELOB
 *
 */

double estimate::compute_mlset_likelihood(multiset* set, cltm_model* model, 
									double *eta_l_v,
									double *var_a_l_v,
									double *var_b_l_v)
{
    double likelihood = 0;
	
    int  v;
	
	for (v = 0; v < model->num_terms; v++){
		likelihood += (log(model->a) - (model->a) * var_a_l_v[v] * var_b_l_v[v]);
		likelihood += (-0.5*log(2*pi)+(0.5-var_a_l_v[v])*(digamma(var_a_l_v[v])+log(var_b_l_v[v]))-0.5 * pow(eta_l_v[v],4));  //*(1/(var_b_l_v[v]*(var_a_l_v[v]-1))));
		likelihood += (var_a_l_v[v]+lgamma(var_a_l_v[v])+var_a_l_v[v]*log(var_b_l_v[v]));
	}
	//printf("ms likelihood: %8.5f \n", likelihood);
    return(likelihood);
}

/*
 * compute global-ELOB
 *
 */

double estimate::compute_global_likelihood(corpus* cps, cltm_model* model, 
									double **var_alpha_phi_k_v,
									double **digamma_subtract_phi_k_v)
{
    double likelihood = 0, var_sum = 0;
    int   k, v;

	likelihood = model->num_topics * lgamma(cps->num_terms * model->alpha_phi)
				- model->num_topics *  cps->num_terms* lgamma(model->alpha_phi);
		
    for (k = 0; k < model->num_topics; k++){
		var_sum=0.0;
		for(v = 0; v < cps->num_terms; v++){
			var_sum += var_alpha_phi_k_v[k][v];
			likelihood += (model->alpha_phi - var_alpha_phi_k_v[k][v])*digamma_subtract_phi_k_v[k][v];
			likelihood += lgamma(var_alpha_phi_k_v[k][v]);
		}
		likelihood = likelihood - lgamma(var_sum);
	}
	//printf("Global likelihood: %8.5f \n", likelihood);
    return(likelihood);
}

double estimate::doc_variational_inference(document* doc, cltm_model* model, 
									double* var_alpha_lamda_lm_i, 
									double* var_alpha_gamma_lm_k,
									double** var_beta_lmn_k,
									double** var_lamda_lmn_i,
									double **digamma_subtract_phi_k_v,
									double *dg_vob_v,
									double log_sum_exp_dg_vob) 
{
    double converged = 1;
    double betasum = 0, lamdasum=0, likelihood = 0;
    double likelihood_old = 0;
	double sum=0.0;
    int  k, i, n, var_iter;
    double* digamma_gam;
	double* digamma_lam;
	double* digamma_subtract_gam;
	

	
	digamma_subtract_gam=(double*)malloc(sizeof(double)*(model->num_topics));
	digamma_gam=(double*)malloc(sizeof(double)*(model->num_topics));
	digamma_lam=(double*)malloc(sizeof(double)*(model->num_lamda));
	
	
    // compute posterior dirichlet

	//initial variational parameters
	// not remove to global
	
	for(k = 0; k < model->num_topics; k++){
		var_alpha_gamma_lm_k[k]= model->alpha_gamma+(doc->total/(((double) model->num_lamda)*((double) model->num_topics)));
		digamma_gam[k]=digamma(var_alpha_gamma_lm_k[k]);
	}
	for(i = 0; i < model->num_lamda; i++){
		var_alpha_lamda_lm_i[i]= model->alpha_lamda+(doc->total/((double) model->num_lamda));
		digamma_lam[i]=digamma(var_alpha_lamda_lm_i[i]);
	}

	
	for (n = 0; n < doc->length; n++){
		for (k = 0; k < model->num_topics; k++){
			var_beta_lmn_k[n][k]=1.0/model->num_topics;
		}
		for(i = 0; i < model->num_lamda; i++){
			var_lamda_lmn_i[n][i]=1.0/model->num_lamda;
		}
	}
	
    // update variationals
    var_iter = 0;

    while (((converged < 0) ||(converged > VAR_DOC_CONVERGED)) && ((var_iter < VAR_MAX_ITER) || (VAR_MAX_ITER == -1)))
    {
		var_iter++;
		likelihood = 0;
		
		for(k = 0; k < model->num_topics; k++){
			digamma_subtract_gam[k]=digamma_subtract(var_alpha_gamma_lm_k, k, model->num_topics);
		}

		for (n = 0; n < doc->length; n++){

			//update var_beta_lmn_k in log space

			betasum = 0;
            for (k = 0; k < model->num_topics; k++){
				//oldbeta[k] = var_beta_lmn_k[n][k];
				/*
				new_var_beta_lmn_k[k]=digamma_subtract_phi_k_v[k][doc->words[n]];
				new_var_beta_lmn_k[k]+=var_lamda_lmn_i[n][1]*digamma_gam[k];
				new_var_beta_lmn_k[k]+=var_lamda_lmn_i[n][0]*dg_topic_k[k];
				*/
				var_beta_lmn_k[n][k] = var_lamda_lmn_i[n][1]*digamma_subtract_phi_k_v[k][doc->words[n]]+digamma_gam[k];
				
                if (k > 0)
                    betasum = log_sum(betasum, var_beta_lmn_k[n][k]);
                else
                    betasum = var_beta_lmn_k[n][k]; // note, phi is in log space
            }

            for (k = 0; k < model->num_topics; k++)
            {
                var_beta_lmn_k[n][k] = exp(var_beta_lmn_k[n][k] - betasum);
            }
			
			//update var_lamda_lmn_i in log space

			lamdasum = 0;
			for(i = 0; i < model->num_lamda; i++){
				//oldlamda[i] = var_lamda_lmn_i[n][i];
				if(i==0){
					var_lamda_lmn_i[n][i]=digamma_lam[i]+dg_vob_v[doc->words[n]]-log_sum_exp_dg_vob;
				}
				else{
					sum=0.0;
					for(k = 0; k < model->num_topics; k++){
						sum += var_beta_lmn_k[n][k]*digamma_subtract_phi_k_v[k][doc->words[n]];
					}
					var_lamda_lmn_i[n][i]=digamma_lam[i]+sum;
				}

                if (i > 0)
                    lamdasum = log_sum(lamdasum, var_lamda_lmn_i[n][i]);
                else
                    lamdasum = var_lamda_lmn_i[n][i]; // note, phi is in log space
			}

			for(i = 0; i < model->num_lamda; i++){
				var_lamda_lmn_i[n][i] = exp(var_lamda_lmn_i[n][i] - lamdasum);
			}

			//update beta
			/*
			for (k = 0; k < model->num_topics; k++)
            {
                var_beta_lmn_k[n][k] = new_var_beta_lmn_k[k];
            }
			*/
		} // end for n : beta and lamda

		
		// update alpha_gamma and alpha_lamda

		for(k = 0; k < model->num_topics; k++){
			var_alpha_gamma_lm_k[k] = model->alpha_gamma;
			for (n = 0; n < doc->length; n++){
				var_alpha_gamma_lm_k[k] += doc->counts[n]* var_beta_lmn_k[n][k];
			}
			digamma_gam[k]=digamma(var_alpha_gamma_lm_k[k]);
		}
		for(i = 0; i < model->num_lamda; i++){
			var_alpha_lamda_lm_i[i] = model->alpha_lamda;
			for (n = 0; n < doc->length; n++){
				var_alpha_lamda_lm_i[i] += doc->counts[n]* var_lamda_lmn_i[n][i];
			}
			digamma_lam[i]=digamma(var_alpha_lamda_lm_i[i]);
		} // end of updates

        likelihood = compute_doc_likelihood( doc, model,
										 var_alpha_lamda_lm_i, 
									     var_alpha_gamma_lm_k,
										 var_beta_lmn_k,
										 var_lamda_lmn_i,
									     digamma_subtract_phi_k_v,
										 dg_vob_v,
									     log_sum_exp_dg_vob);
        assert(!isnan(likelihood));
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;

        //printf("doc %s elob: %8.5f %1.3e\n",doc->docid.c_str(), likelihood, converged);

    } //end while

	//printf("Final doc %s elob: %8.5f %1.3e\n", doc->docid.c_str(), likelihood, converged);
	//free malloc
	
	free (digamma_gam);
	free (digamma_lam);
	free (digamma_subtract_gam);
	//free (new_var_beta_lmn_k);
    return likelihood;
}


double estimate::mlset_variational_inference(multiset* set, cltm_model* model, 
									double** var_alpha_lamda_lm_i, 
									double** var_alpha_gamma_lm_k,
									double*** var_beta_lmn_k,
									double*** var_lamda_lmn_i,
									double **digamma_subtract_phi_k_v,
									double *var_theta_v,
									double *eta_l_v,
									double *var_a_l_v,
									double *var_b_l_v)
{
    double converged = 1;
    double likelihood = 0;
    double likelihood_old = 0;
	double sum=0.0,sum1=0.0;
	double log_sum_exp_dg_vob=0.0;
    int  m, v, var_iter;
	
	double* dg_vob_v;
	//double* new_eta_l_k;
	//double *new_var_a_l_k;
	//double *new_var_b_l_k;

    // compute posterior dirichlet

	//prepare some fixed cals from the global variations
	dg_vob_v=(double*)malloc(sizeof(double)*(model->num_terms));
	//new_eta_l_k=(double*)malloc(sizeof(double)*(model->num_topics));
	//new_var_a_l_k=(double*)malloc(sizeof(double)*(model->num_topics));
	//new_var_b_l_k=(double*)malloc(sizeof(double)*(model->num_topics));


	
    // update variationals
    var_iter = 0;


    while (( fabs(converged) > VAR_SET_CONVERGED) &&
           ((var_iter < VAR_SET_MAX_ITER) || (VAR_SET_MAX_ITER == -1)))
    {
		var_iter++;
		likelihood = 0;

		

		log_sum_exp_dg_vob=0.0;
		for(v = 0; v < model->num_terms; v++){
			dg_vob_v[v]=var_theta_v[v]+eta_l_v[v];
			log_sum_exp_dg_vob+=exp(dg_vob_v[v]);
		}
		log_sum_exp_dg_vob=log(log_sum_exp_dg_vob);


		// need to update first, contain the initialization for beta and lamda
		// update local doc variables
		for (m = 0; m < set->size; m++){
				likelihood +=doc_variational_inference(&(set->docs[m]), model, 
										var_alpha_lamda_lm_i[m], 
										var_alpha_gamma_lm_k[m],
										var_beta_lmn_k[m],
										var_lamda_lmn_i[m],
										digamma_subtract_phi_k_v,
										dg_vob_v,
										log_sum_exp_dg_vob);

		}
		
        likelihood += compute_mlset_likelihood(set, model, 
									eta_l_v,
									var_a_l_v,
									var_b_l_v);
		
		////////////////////////
		//update a, b, eta
		///
		map<char,void*> paras_map;
		paras_map['a']=var_theta_v;
		paras_map['b']=var_lamda_lmn_i;
		//paras_map['c']=var_beta_lmn_k;
		paras_map['d']=set;

		eta_l_v=opt_eta_gsl(eta_l_v,&paras_map,model->num_terms);
		///
		for(v = 0; v < model->num_terms; v++){
			var_a_l_v[v]=opt_a(var_a_l_v[v],var_b_l_v[v],eta_l_v[v]); 
			var_b_l_v[v]=pow(eta_l_v[v],2)/(var_a_l_v[v]-1);
		}
		
		/* sth wrong 
		map<char,void*> paras_mapa;
		paras_mapa['a']=eta_l_k;
		paras_mapa['b']=var_b_l_k;
		var_a_l_k=opt_a_gsl(var_a_l_k,&paras_mapa,model->num_topics);

		map<char,void*> paras_mapb;
		paras_mapb['a']=eta_l_k;
		paras_mapb['b']=var_a_l_k;
		var_b_l_k=opt_b_gsl(var_b_l_k,&paras_mapb,model->num_topics);
		*/
		//end of update a, b, eta
		////////////////////////
						
        assert(!isnan(likelihood));
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;

        // printf("[LDA INF] %8.5f %1.3e\n", likelihood, converged);
		//printf("Set elob: %8.5f %1.3e\n", likelihood, converged);
    } //end while

	printf("Final set elob: %8.5f %1.3e\n", likelihood, converged);

	//free
	//free (new_eta_l_k);
	//free (new_var_a_l_k);
	//free (new_var_b_l_k);
	free (dg_vob_v);

	return(likelihood);
}

//Problem: global elob decrease every iteration

void balance_theta_eta(double *var_theta_v,double **eta_l_v,int num_terms, int  num_multisets){
	int i, v,l,min_index;
	for(v = 0; v < num_terms; v++){
		min_index=0;
		for(l = 0; l < num_multisets; l++){
			if(eta_l_v[l][v]<=0)
				break;
			if(eta_l_v[l][v]<eta_l_v[min_index][v]){
				min_index=l;
			}
		}
		if(l==num_multisets){
			var_theta_v[v]+=eta_l_v[min_index][v];
			for(i = 0;i < num_multisets; i++){
				eta_l_v[i][v]-=eta_l_v[min_index][v];
			}
		}
	}
}
double estimate::global_variational_inference(corpus* cps, cltm_model* model, 
									double*** var_alpha_lamda_lm_i, 
									double*** var_alpha_gamma_lm_k,
									double**** var_beta_lmn_k,
									double**** var_lamda_lmn_i,
									double **var_alpha_phi_k_v,
									double *var_theta_v,
									double **eta_l_v,
									double **var_a_l_v,
									double **var_b_l_v,
									char* directory)
{
	
	double converged = 1;
    double likelihood = 0;
    double likelihood_old = 0;
	double** digamma_subtract_phi_k_v;
	
    int  l,  m, v, k,  n,  iter;
	char filename[100];
	
    // compute posterior dirichlet

	// initial variational parameters
	
	for(k = 0; k < model->num_topics; k++){
		for(v = 0; v < cps->num_terms; v++){
			var_alpha_phi_k_v[k][v]= model->alpha_phi + (cps->term_counts[v]/((double) model->num_topics));
		}
	}

	srand(time(0));

	for(v = 0; v < cps->num_terms; v++){
		var_theta_v[v]=((rand()%10)+1)/10.0;
	}
	

	//allocate digamma_subtract;
	digamma_subtract_phi_k_v=(double**)malloc(sizeof(double*)*(model->num_topics));
	for(k = 0; k < model->num_topics; k++){
		digamma_subtract_phi_k_v[k]=(double*)malloc(sizeof(double)*(model->num_terms));
		for(v = 0; v < model->num_terms; v++){
			digamma_subtract_phi_k_v[k][v]=digamma_subtract(var_alpha_phi_k_v[k],v,model->num_terms);
		}
	}
	std::default_random_engine generator;
	for (l = 0; l < cps->num_multisets; l++){
		for(v = 0; v < cps->num_terms; v++){
			var_a_l_v[l][v]=10.0;
			var_b_l_v[l][v]=5.0;
			std::gamma_distribution<double> gamma_dir(var_a_l_v[l][v],var_b_l_v[l][v]);
			double variance = gamma_dir(generator);
			std::normal_distribution<double> normal_dir(0.0,variance);
			//eta_l_k[l][k] = normal_dir(generator); //global converged, word prob too small
			eta_l_v[l][v] = 10.0; //local converged, word prob ok
			//eta_l_k[l][k] = rand()%10;  //local global not converged
			//eta_l_k[l][k] = normal_dir(generator);//local global not converged, word prob too small
		}
	}
	/*
	for(l = 0; l < cps->num_multisets; l++){
		for(m = 0; m < cps->mlsets[l].size; m++){
			for(k = 0; k < model->num_topics; k++){
				var_alpha_gamma_lm_k[l][m][k]= model->alpha_gamma+(cps->mlsets[l].docs[m].total/(((double) model->num_lamda)*((double) model->num_topics)));
		
			}
			for(i = 0; i < model->num_lamda; i++){
				var_alpha_lamda_lm_i[l][m][i]= model->alpha_lamda+(cps->mlsets[l].docs[m].total/((double) model->num_lamda));
		
			}
			for(n=0;n< cps->mlsets[l].docs[m].length;n++){
				for (k = 0; k < model->num_topics; k++){
					var_beta_lmn_k[l][m][n][k]=1.0/model->num_topics;
				}
				for(i = 0; i < model->num_lamda; i++){
					var_lamda_lmn_i[l][m][n][i]=1.0/model->num_lamda;
				}
			}
		}
	}
	*/

    // compute posterior dirichlet

    iter = 0;
	while (((converged < 0) || (converged > EM_CONVERGED)) && ((iter < EM_MAX_ITER) || (VAR_MAX_ITER == -1)))
    // while (((converged < 0) || (converged > EM_CONVERGED) || (iter <= 2)) && (iter <= EM_MAX_ITER))
    {
		iter++;
		printf("**** em iteration %d ****\n", iter);
		likelihood = 0;
		

		for (l = 0; l < cps->num_multisets; l++){
			likelihood +=mlset_variational_inference(&(cps->mlsets[l]), model, 
										var_alpha_lamda_lm_i[l], 
										var_alpha_gamma_lm_k[l],
										var_beta_lmn_k[l],
										var_lamda_lmn_i[l],
										digamma_subtract_phi_k_v,
										var_theta_v,
										eta_l_v[l],
										var_a_l_v[l],
										var_b_l_v[l]);

		}


		// compute elob
		
       likelihood += compute_global_likelihood(cps, model, 
									var_alpha_phi_k_v,
									digamma_subtract_phi_k_v);
		
	   
	   //update alpha_phi, alpha_theta
		for(k = 0; k < model->num_topics; k++){
			for(v = 0; v < cps->num_terms; v++){
				var_alpha_phi_k_v[k][v]= model->alpha_phi;
				//var_theta_k[k]=model->alpha_theta;
			}
		}

		l=0;m=0;n=0;
		cout<<"theta: ";
		for(v = 0; v < 10; v++){
			cout<<var_theta_v[v]<<"  ";
		}
		cout<<endl;
		cout<<"eta: ";
		for(v = 0; v < 10; v++){
			cout<<eta_l_v[l][v]<<"  ";
		}
		cout<<endl;
		cout<<"beta: ";
		for(k = 0; k < 10; k++){
			cout<<var_beta_lmn_k[l][m][n][k]<<"  ";
		}
		cout<<endl;
	
		for(l = 0; l < cps->num_multisets; l++){
		
			for(m = 0; m < (cps->mlsets[l]).size; m++){

				for(n = 0; n <  (cps->mlsets[l]).docs[m].length; n++){
						
					for(k = 0; k < model->num_topics; k++){
						//var_theta_k[k] += (cps->mlsets[l].docs[m].counts[n]) * var_beta_lmn_k[l][m][n][k] * var_lamda_lmn_i[l][m][n][0];
						var_alpha_phi_k_v[k][cps->mlsets[l].docs[m].words[n]] += (cps->mlsets[l].docs[m].counts[n] *var_lamda_lmn_i[l][m][n][1] *var_beta_lmn_k[l][m][n][k]);

					}
						
				}
			}
		}
		// compute digamma_subtract
		
		for(k = 0; k < model->num_topics; k++){
			for(v = 0; v < model->num_terms; v++){
				digamma_subtract_phi_k_v[k][v]=digamma_subtract(var_alpha_phi_k_v[k],v,model->num_terms);
			}
		}


		/*
		map<char,void*> paras_map;
		paras_map['a']=eta_l_k;
		paras_map['b']=var_lamda_lmn_i;
		paras_map['c']=var_beta_lmn_k;
		paras_map['d']=cps;

		var_theta_k=opt_theta_gsl(var_theta_k,&paras_map,model->num_topics);
		*/

		//end of update alpha_phi, alpha_theta

		balance_theta_eta(var_theta_v,eta_l_v,cps->num_terms, cps->num_multisets);
		 
		///////////////////////////////////		

        assert(!isnan(likelihood));
        converged = (likelihood_old - likelihood) / likelihood_old;
        likelihood_old = likelihood;

		 printf("Global elob: %8.5f %1.3e\n", likelihood, converged);
        // printf("[LDA INF] %8.5f %1.3e\n", likelihood, converged);

		  if ((iter % LAG) == 0)
		  {
            sprintf(filename,"%s\\%03d.theta0",directory, iter);
            save_theta0(filename, var_theta_v, cps->num_terms);
			sprintf(filename,"%s\\%03d.thetag",directory, iter);
			save_thetag(filename, var_theta_v,eta_l_v, cps->num_multisets,cps->num_terms);
			sprintf(filename,"%s\\%03d.alphi",directory, iter);
			save_alpha_phi(filename, var_alpha_phi_k_v, model->num_topics, cps->num_terms);
        }
	


    }//end while

	 printf("Final global elob: %8.5f %1.3e\n", likelihood, converged);
	//free malloc
	for(k = 0; k < model->num_topics; k++){
		free (digamma_subtract_phi_k_v[k]);
	}
	free (digamma_subtract_phi_k_v);

	
	//free (T_l);

	//free (T1_l);

    return(likelihood);
}


void estimate::run_vi(corpus* corpus,char* directory)
{

    int l, m, n, k;
    cltm_model *model = NULL;
	double **var_alpha_phi_k_v;
	double ***var_alpha_lamda_lm_i;
	double ***var_alpha_gamma_lm_k;
	double *var_theta_v;
	double ****var_beta_lmn_k;
	double ****var_lamda_lmn_i;
	double **var_a_l_v;
	double **var_b_l_v;
	double **eta_l_v;
	char filename[100];

	// initialize model

	model=(cltm_model*)malloc(sizeof(cltm_model));

	//double phi,double lamda,double gamma,double theta,int nterm,int ntopics
	model->model_init(0.01,1.0,0.1,1.0,corpus->num_terms,10);

    // allocate variational parameters
	
	var_alpha_phi_k_v = (double**)malloc(sizeof(double*)*(model->num_topics));
    for(k = 0; k < model->num_topics; k++){
		var_alpha_phi_k_v[k]=(double*)malloc(sizeof(double) * (corpus->num_terms));
	}
    
	var_theta_v = (double*)malloc(sizeof(double)*(corpus->num_terms));

	var_alpha_lamda_lm_i = (double***)malloc(sizeof(double**)*(corpus->num_multisets));
    for(l = 0; l < corpus->num_multisets; l++){
		var_alpha_lamda_lm_i[l]=(double**)malloc(sizeof(double*) * (corpus->mlsets[l].size));
		for(m = 0; m < corpus->mlsets[l].size; m++){
			var_alpha_lamda_lm_i[l][m]=(double*)malloc(sizeof(double) * (model->num_lamda));
		}
	}

	var_alpha_gamma_lm_k = (double***)malloc(sizeof(double**)*(corpus->num_multisets));
    for(l = 0; l < corpus->num_multisets; l++){
		var_alpha_gamma_lm_k[l]=(double**)malloc(sizeof(double*) * (corpus->mlsets[l].size));
		for(m = 0; m < corpus->mlsets[l].size; m++){
			var_alpha_gamma_lm_k[l][m]=(double*)malloc(sizeof(double) * (model->num_topics));
		}
	}

	var_beta_lmn_k = (double****)malloc(sizeof(double***)*(corpus->num_multisets));
    for(l = 0; l < corpus->num_multisets; l++){
		var_beta_lmn_k[l]=(double***)malloc(sizeof(double**) * (corpus->mlsets[l].size));
		for(m = 0; m < corpus->mlsets[l].size; m++){
			var_beta_lmn_k[l][m]=(double**)malloc(sizeof(double*) * (corpus->mlsets[l].docs[m].length));
			for(n=0;n<corpus->mlsets[l].docs[m].length;n++){
				var_beta_lmn_k[l][m][n]=(double*)malloc(sizeof(double) * (model->num_topics));
			}
		}
	}

	var_lamda_lmn_i = (double****)malloc(sizeof(double***)*(corpus->num_multisets));
    for(l = 0; l < corpus->num_multisets; l++){
		var_lamda_lmn_i[l]=(double***)malloc(sizeof(double**) * (corpus->mlsets[l].size));
		for(m = 0; m < corpus->mlsets[l].size; m++){
			var_lamda_lmn_i[l][m]=(double**)malloc(sizeof(double*) * (corpus->mlsets[l].docs[m].length));
			for(n=0;n<corpus->mlsets[l].docs[m].length;n++){
				var_lamda_lmn_i[l][m][n]=(double*)malloc(sizeof(double) * (model->num_lamda));
			}
		}
	}

	var_a_l_v = (double**)malloc(sizeof(double*)*(corpus->num_multisets));
    for(l = 0; l < corpus->num_multisets; l++){
		var_a_l_v[l]=(double*)malloc(sizeof(double) * (corpus->num_terms));
	}
	var_b_l_v = (double**)malloc(sizeof(double*)*(corpus->num_multisets));
    for(l = 0; l < corpus->num_multisets; l++){
		var_b_l_v[l]=(double*)malloc(sizeof(double) * (corpus->num_terms));
	}
	eta_l_v = (double**)malloc(sizeof(double*)*(corpus->num_multisets));
    for(l = 0; l < corpus->num_multisets; l++){
		eta_l_v[l]=(double*)malloc(sizeof(double) * (corpus->num_terms));
	}

    // run vi

	double elob=global_variational_inference(corpus,  model, 
									var_alpha_lamda_lm_i, 
									var_alpha_gamma_lm_k,
									var_beta_lmn_k,
									var_lamda_lmn_i,
									var_alpha_phi_k_v,
									var_theta_v,
									eta_l_v,
									var_a_l_v,
									var_b_l_v,directory);

	cout<<"Final ELOB: "<<elob<<endl;

	sprintf(filename,"%s\\final.theta0",directory);
    save_theta0(filename, var_theta_v,corpus->num_terms);
	sprintf(filename,"%s\\final.thetag",directory);
	save_thetag(filename, var_theta_v,eta_l_v, corpus->num_multisets, corpus->num_terms);
	sprintf(filename,"%s\\final.alphi",directory);
	save_alpha_phi(filename, var_alpha_phi_k_v, model->num_topics, corpus->num_terms);

	//free
	 
    for(k = 0; k < model->num_topics; k++){
		free (var_alpha_phi_k_v[k]);
	}
    free (var_alpha_phi_k_v);
	free (var_theta_v);

    for(l = 0; l < corpus->num_multisets; l++){
		
		for(m = 0; m < corpus->mlsets[l].size; m++){
			free (var_alpha_lamda_lm_i[l][m]);
		}
		free (var_alpha_lamda_lm_i[l]);
	}
	free (var_alpha_lamda_lm_i);

    for(l = 0; l < corpus->num_multisets; l++){
		for(m = 0; m < corpus->mlsets[l].size; m++){
			free (var_alpha_gamma_lm_k[l][m]);
		}
		free (var_alpha_gamma_lm_k[l]);
	}
	free (var_alpha_gamma_lm_k);

	
    for(l = 0; l < corpus->num_multisets; l++){		
		for(m = 0; m < corpus->mlsets[l].size; m++){			
			for(n=0;n<corpus->mlsets[l].docs[m].length;n++){
				 free (var_beta_lmn_k[l][m][n]);
			}
			free (var_beta_lmn_k[l][m]);
		}
		free (var_beta_lmn_k[l]);
	}
	free (var_beta_lmn_k);

	
    for(l = 0; l < corpus->num_multisets; l++){
		for(m = 0; m < corpus->mlsets[l].size; m++){
			for(n=0;n<corpus->mlsets[l].docs[m].length;n++){
				free (var_lamda_lmn_i[l][m][n]);
			}
			free (var_lamda_lmn_i[l][m]);
		}
		free (var_lamda_lmn_i[l]);
	}
	free (var_lamda_lmn_i);

    for(l = 0; l < corpus->num_multisets; l++){
		free (var_a_l_v[l]);
		free (var_b_l_v[l]);
		free (eta_l_v[l]);
	}
	free (var_a_l_v);
	free (var_b_l_v);
	free (eta_l_v);

}
