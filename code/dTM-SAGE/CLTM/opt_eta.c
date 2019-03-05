// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)

// This file is part of LDA-C.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#include "opt_eta.h"


/*
 * objective function and its derivatives
 *
 */

double d_eta(double eta, double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double a,double b, double T,
			   multiset* set,int cur_k)
{ 
	int m,n;

	double res=0.0;
	for(m = 0; m < set->size; m++){
		for(n=0;n<set->docs[m].length;n++){
			res+=(var_lamda_lmn_i[m][n][0]*var_beta_lmn_k[m][n][cur_k]);
		}
	}
	return (1-T)*res+(-eta*(1/(b*(a-1))));
}

double d2_eta(double eta, double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double a,double b, double T,
			   multiset* set,int cur_k)
{   int m,n;
	double res=0.0;
	for(m = 0; m < set->size; m++){
		for(n=0;n<set->docs[m].length;n++){
			res+=var_lamda_lmn_i[m][n][0]*var_beta_lmn_k[m][n][cur_k];
		}
	}
	return res*(pow(T,2)-T)+(-1/(b*(a-1)));
}

/*
 * newtons method
 *
 */
/*
									eta_l_k[k],var_beta_lmn_k,
									var_lamda_lmn_i,
									var_a_l_k[k],
									var_b_l_k[k],
									T,set,model,k
									*/

double opt_eta(double eta_init,double*** var_beta_lmn_k,
			   double*** var_lamda_lmn_i,
			   double a,double b, double T,
			   multiset* set,int cur_k)
{
    double eta, log_eta;
    double df, d2f;
    int iter = 0;

    log_eta = log(eta_init);
    do
    {
        iter++;
        eta = exp(log_eta);
        if (isnan(eta))
        {
            eta_init = eta_init * 10;
            printf("warning : alpha is nan; new init = %5.5f\n", eta_init);
            eta =eta_init;
            log_eta = log(eta);
        }
   
        df = d_eta(eta, var_beta_lmn_k,
			   var_lamda_lmn_i,
			   a,b,T,set,cur_k);
        d2f = d2_eta( eta,  var_beta_lmn_k,
			    var_lamda_lmn_i,
			   a,b,T,set,cur_k);

        log_eta = log_eta - df/(d2f * eta + df);
       // printf("alpha maximization : %5.5f   %5.5f\n", f, df);
    }
    while ((fabs(df) > NEWTON_THRESH) && (iter < MAX_ETA_ITER));
    return(exp(log_eta));
}
