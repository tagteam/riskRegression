#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>       /* isnan */

#ifndef ARMA_WRAP_H
#define ARMA_WRAP_H

#define ARMA_NO_DEBUG

#define GET_VARIABLE_NAME(Variable) (#Variable)

#ifndef ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_OPENMP 1
#endif

#include <RcppArmadillo.h>

void checkNAs(Rcpp::NumericVector& vec, std::string var_name);
void checkNAs(double val, std::string var_name);
void compareLengths(Rcpp::NumericVector& vec1, Rcpp::NumericVector& vec2);

#endif
