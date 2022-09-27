#include "arma-wrap.h"

using namespace arma;

void checkNAs(Rcpp::NumericVector& vec, std::string var_name){
  for (int i=0; i< vec.size(); i++) {
    if (R_IsNA(vec[i])) {
      Rcpp::stop("Missing values in variable %i. ", var_name); // Calls R function for stopping
    }
  }
}

void checkNAs(double val, std::string var_name){
  if (isnan(val)){
    Rcpp::stop("Missing values in variable %i. ", var_name);
  }
}

void compareLengths(Rcpp::NumericVector& vec1, Rcpp::NumericVector& vec2){
  if (vec1.size() != vec2.size()){
    Rcpp::stop("Some vectors have unequal lengths. ");
  }
}