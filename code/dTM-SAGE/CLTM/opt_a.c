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

#include "opt_a.h"

/*
 * objective function and its derivatives
 *
 */


double dl_a(double a, double b, double eta)
{ return (0.5+a)*trigamma(a) - 0.5*pow(eta,2)*(1/b)*(1/pow((a-1),2)) - 1;}

double d2l_a(double a, double b, double eta)
{ return (0.5+a)*quadgamma(a)- pow(eta,2)*(1/b)*(1/pow((a-1),3));}


/*
 * newtons method
 *
 */

double opt_a(double init_a, double b, double eta)
{
    double a, log_a;
    double df, d2f;
    int iter = 0;

    log_a = log(init_a);
    do
    {
        iter++;
        a = exp(log_a);
        if (isnan(a))
        {
            init_a = init_a * 10;
            printf("warning : alpha is nan; new init = %5.5f\n", init_a);
            a = init_a;
            log_a = log(a);
        }
        df = dl_a(a, b, eta);
        d2f = d2l_a(a, b, eta);
        log_a = log_a - df/(d2f * a + df);
        printf("alpha maximization : %5.5f   %5.5f\n", f, df);
    }
    while ((fabs(df) > NEWTON_THRESH) && (iter < MAX_A_ITER));
    return(exp(log_a));
}
