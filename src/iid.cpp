// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

double calcS0_cpp(double t, int n, const NumericVector& eventtime, const NumericVector& eXb);
NumericVector calcS1_cpp(double t, int n, int p, const NumericVector& eventtime, const NumericVector& eXb, const arma::mat& X);

double calcS0_cpp(double t, int n,
                  const NumericVector& eventtime,
                  const NumericVector& eXb){
  
  double S0=0;
  
  for(int iterObs=0;iterObs<n;iterObs++){
    if(t<=eventtime[iterObs]){
      S0 += eXb[iterObs];
    }
  }
  
  return(S0);
  
}

NumericVector calcS1_cpp(double t, int n, int p,
                         const NumericVector& eventtime,
                         const NumericVector& eXb,
                         const arma::mat& X){
  
  NumericVector S1(p,0.0);
  for(int iterObs=0;iterObs<n;iterObs++){
    for(int iterX=0;iterX<p;iterX++){
      if(t<=eventtime[iterObs]){
        S1[iterX] += X(iterObs,iterX) * eXb[iterObs];
      }
    }
  }
  
  return(S1);
  
}

// [[Rcpp::export]]
List calcEstrata_cpp(double tau,
                     const IntegerVector& strata,
                     const IntegerVector& status,
                     int n, int p, int nStrata,
                     const NumericVector& eventtime,
                     const NumericVector& eXb,
                     const arma::mat& X){
  
  arma::colvec strata2vec = as<colvec>(strata);
  
  // store results
  vector<vec> lsS0(nStrata);
  vector<arma::mat> lsS1(nStrata);
  vector<arma::mat> lsE(nStrata);
  vector<NumericVector> lsUtime1(nStrata);
  
  // tempo
  NumericVector sub_eventtime, ssub_eventtime, sub_statut, sub_eXb;
  arma::mat sub_X;
  int n_iterUtime1;
  
  double S0;
  NumericVector S1, E;
    
  for(int iterStrata=0; iterStrata<nStrata; iterStrata++){
    R_CheckUserInterrupt();
    
    Rcout << iterStrata << endl;
    // subset input according to strata
    sub_eventtime = eventtime[strata==iterStrata];
    sub_statut = status[strata==iterStrata];
    sub_eXb = eXb[strata==iterStrata];
    sub_X = X.rows(find(strata2vec==iterStrata));
    
    Rcout << "a" << endl;
    // vector of unique time for events
    ssub_eventtime = sub_eventtime[sub_statut>0];
    lsUtime1[iterStrata] = unique(ssub_eventtime);
    n_iterUtime1 = lsUtime1[iterStrata].size();
    
    Rcout << "b: " <<  n_iterUtime1 << " " << p << endl;
    // prepare output
    // lsS0[1].resize(n_iterUtime1);
    // lsS1[1].resize(n_iterUtime1, p);
    // lsE[1].resize(n_iterUtime1, p);

    Rcout << "c" << endl;
    // for(int iterUTime1=0; iterUTime1<n_iterUtime1; iterUTime1++){
    //   if(lsUtime1[iterStrata][iterUTime1]>tau){
    //     
    //     lsS0[iterUTime1][iterUTime1] = 0;
    //     for(int iterP=0; iterP<p; iterP++){
    //       lsS1[iterUTime1](iterUTime1,iterP) = 0;
    //       lsE[iterUTime1](iterUTime1,iterP) = 0;
    //     }
    //     
    //   }else{
    //     
    //     S0 = calcS0_cpp(lsUtime1[iterStrata][iterUTime1], n, sub_eventtime, sub_eXb);
    //     S1 = calcS1_cpp(lsUtime1[iterStrata][iterUTime1], n, p, sub_eventtime, sub_eXb, sub_X);
    //     E = S1/S0;
    //     if(S0==0){std::fill(E.begin(), E.end(), 0.0);}
    //     lsS0[iterUTime1][iterUTime1] = S0;
    //     lsS1[iterUTime1].row(iterUTime1) = as<rowvec>(S1);
    //     lsE[iterUTime1].row(iterUTime1) = as<rowvec>(E);
    //     
    //   }
    // }
    
  }
  
  return(List::create(Named("Utime1")  = lsUtime1,
                      Named("S0")  = lsS0,
                      Named("S1")  = lsS1,
                      Named("E")  = lsE));
}

// [[Rcpp::export]]
arma::mat calcU_cpp(const arma::mat& newX, const NumericVector& newStatus, int newN,
                    const IntegerVector& IndexNewT, const arma::mat& ENewT, 
                    int p, bool aggregate){
  
  // NumericVector U1_time = newT[newStatus>0];//unique();
  // U1_time = unique(U1_time);
  // IntegerVector n2u1 = match(newT, U1_time)-1;
  
  // compute the score
  arma::mat score(newN, p);
  score.fill(0.0);
  
  for(int iterObs=0;iterObs<newN;iterObs++){
    if(newStatus[iterObs]==0){continue;}
    score.row(iterObs) = newX.row(iterObs) - ENewT.row( IndexNewT[iterObs] );
  }
  
  if(aggregate){
    score <- sum(score,0);
  }
  
  // export
  return(score);
}
