#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef ARMA_WRAP_H
#define ARMA_WRAP_H

#define ARMA_NO_DEBUG

#ifndef ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_OPENMP 1
#endif

#include <RcppArmadillo.h>

#endif
