#include <RcppArmadillo.h>

void getInfluenceFunctionKM(Rcpp::NumericVector& time, Rcpp::NumericVector& status,arma::vec& atrisk,arma::vec& MC_term2,arma::uvec& sindex,arma::vec& utime);

