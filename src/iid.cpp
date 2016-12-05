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
List calcE_cpp(double t, int n, int p,
               const NumericVector& eventtime,
               const NumericVector& eXb,
               const arma::mat& X){
  
  double S0 = calcS0_cpp(t, n, eventtime, eXb);
  NumericVector S1 = calcS1_cpp(t, n, p, 
                                eventtime, eXb, X);
  NumericVector E = S1/S0;
  if(S0==0){
    std::fill(E.begin(), E.end(), 0.0);
  }
  
  return(List::create(Named("E")  = E,
                      Named("S1")  = S1,
                      Named("S0")  = S0));
}

// // [[Rcpp::export]]
// List calcEstrata_cpp(double tau,
//                      const vector<IntegerVector>& indexStrata,
//                      const IntegerVector& status,
//                      int p, int nStrata,
//                      const NumericVector& eventtime,
//                      const NumericVector& eXb,
//                      const arma::mat& X, bool add0){
// 
//   double valDefault = 0; // NA_REAL
//   double maxTime = eventtime[eventtime.size()-1];
//   
//   // store results
//   vector<vec> lsS0(nStrata);
//   vector<arma::mat> lsS1(nStrata);
//   vector<arma::mat> lsE(nStrata);
//   vector<NumericVector> lsUtime1(nStrata);
//   NumericVector n(nStrata);
//   
//   // tempo
//   NumericVector sub_eventtime, ssub_eventtime, sub_statut, sub_eXb, timeTempo;
//   arma::mat sub_X;
//   int n_iterUtime1;
//   arma::vec strataVec;
//   
//   double S0;
//   NumericVector S1, E;
//    
//   for(int iterStrata=0; iterStrata<nStrata; iterStrata++){
//     R_CheckUserInterrupt();
//     
//     n[iterStrata] = indexStrata[iterStrata].size();
//     
//     // subset input according to strata
//     sub_eventtime = eventtime[indexStrata[iterStrata]];
//     sub_statut = status[indexStrata[iterStrata]];
//     sub_eXb = eXb[indexStrata[iterStrata]];
//     sub_X = X.rows(as<uvec>(indexStrata[iterStrata]));
//     
//     // vector of unique time for events
//     ssub_eventtime = sub_eventtime[sub_statut>0];
//     timeTempo = unique(ssub_eventtime);
//     std::sort(timeTempo.begin(),timeTempo.end());
//     if(add0){timeTempo.push_back(maxTime+1e-12);}
//     n_iterUtime1 = timeTempo.size();
//     
//     // prepare output
//     lsS0[iterStrata].set_size(n_iterUtime1);
//     lsS1[iterStrata].set_size(n_iterUtime1,p);
//     lsE[iterStrata].set_size(n_iterUtime1,p);
// 
//     for(int iterUTime1=0; iterUTime1<n_iterUtime1; iterUTime1++){
// 
//       if(timeTempo[iterUTime1] > tau){ // after time of interest
//         lsS0[iterStrata][iterUTime1] = valDefault;
//         for(int iterP=0; iterP<p; iterP++){
//           lsS1[iterStrata](iterUTime1,iterP) = valDefault;
//           lsE[iterStrata](iterUTime1,iterP) = valDefault;
//         }
// 
//       }else{ // could be improved 
//         S0 = calcS0_cpp(timeTempo[iterUTime1], n[iterStrata], sub_eventtime, sub_eXb);
//         S1 = calcS1_cpp(timeTempo[iterUTime1], n[iterStrata], p, sub_eventtime, sub_eXb, sub_X);
//         E = S1/S0;
//         if(S0==0){std::fill(E.begin(), E.end(), 0.0);}
//         lsS0[iterStrata][iterUTime1] = S0;
//         lsS1[iterStrata].row(iterUTime1) = as<rowvec>(S1);
//         lsE[iterStrata].row(iterUTime1) = as<rowvec>(E);
//       }
//     }
//     
//     lsUtime1[iterStrata] = timeTempo;
//   }
//   
//   return(List::create(Named("Utime1")  = lsUtime1,
//                       Named("S0")  = lsS0,
//                       Named("S1")  = lsS1,
//                       Named("E")  = lsE,
//                       Named("n")  = n));
// }

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
