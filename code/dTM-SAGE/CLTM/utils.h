#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>

//#include <sys/stat.h>
//#include <sys/types.h>
#include <vector>

bool isnan(double x);
double sum_array(double* vec, int n);
double log_sum(double log_a, double log_b);
double trigamma(double x);
double digamma(double x);
double lgamma(double x);
double quadgamma(double x);
int argmax(double* x, int n);

#endif
