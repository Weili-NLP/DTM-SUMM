#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "utils.h"

#define NEWTON_THRESH 1e-5
#define MAX_A_ITER 1000


double dl_a(double a, double b, double eta);
double d2l_a(double a, double b, double eta);
double opt_a(double init_a, double b, double eta);

//void maximize_alpha(double** gamma, cltm_model* model, int num_docs);


