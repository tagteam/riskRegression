// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

double calcS0_cpp(double t, int n, const NumericVector& eventtime, const NumericVector& eXb);
NumericVector calcS1_cpp(double t, int n, int p, const NumericVector& eventtime, const NumericVector& eXb, const arma::mat& X);
List calcE_cpp(double t, int n, int p, const NumericVector& eventtime, const NumericVector& eXb, const arma::mat& X);

// [[Rcpp::export]]
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

// [[Rcpp::export]]
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


// // [[Rcpp::export]]
// arma::mat predictCIF_cpp(const std::vector<arma::mat>& hazard, 
//                          const std::vector<arma::mat>& cumHazard, 
//                          const arma::mat& eXb_h, 
//                          const arma::mat& eXb_cumH, 
//                          const arma::mat& strata, 
//                          const std::vector<double>& newtimes, 
//                          const std::vector<double>& etimes, 
//                          const std::vector<double>& etimeMax, 
//                          double t0,
//                          int nTimes, 
//                          int nNewTimes, 
//                          int nData, 
//                          int cause, 
//                          int nCause){
// 
//   arma::mat pred_CIF(nData, nNewTimes);
//   pred_CIF.fill(0);
//   
//   double hazard_it;
//   double AllcumHazard_it;
//   double survival_t0=1;
//   int iterP; // index of the prediction time
//   rowvec strataI(nCause);
//   
//   for(int iterI=0 ; iterI<nData; iterI++){ // index of the patient
//     R_CheckUserInterrupt();
//     
//     iterP = 0;
//     strataI = strata.row(iterI);
//     
//     for(int iterT=0 ; iterT<nTimes; iterT++){ // index of the time in the integral (event time number)
//       // update position 
//       while(iterP < nNewTimes && newtimes[iterP]<etimes[iterT]){
//         iterP++;
// 	pred_CIF(iterI,iterP) = pred_CIF(iterI,iterP-1);
//       }
//       if(iterP >= nNewTimes){break;}
//       // get hazard for the cause of interest
//       hazard_it = hazard[cause](iterT,strataI[cause])*eXb_h(iterI,cause);
//       
//       // sum all cumHazard for all causes and exp the result
//       AllcumHazard_it = 0; 
//       for(int iterC=0 ; iterC<nCause; iterC++){
//         AllcumHazard_it += cumHazard[iterC](iterT,strataI[iterC])*eXb_cumH(iterI,iterC);
//       }
//       
//       // update the integral
//       if(R_IsNA(t0)){
//         pred_CIF(iterI,iterP) += exp(-AllcumHazard_it) * hazard_it;  
//       }else{// [only for conditional CIF]
//         
//         if(etimes[iterT]<t0 && iterT<(nTimes+1) && etimes[iterT+1]>=t0){
//           survival_t0 = exp(-AllcumHazard_it); // get the survival up to t0
//         }
//         
//         if(etimes[iterT] >= t0){ // not needed  newtimes[iterP]>=t0  because newtimes >= etimes see update position above 
//           // Rcout << hazard_it << " ("<< iterP << ","<< etimes[iterT] << ","<< newtimes[iterP] << ")";
//           pred_CIF(iterI,iterP) += exp(-AllcumHazard_it) * hazard_it / survival_t0;
//         }
//         
//       }
//     }
//     // Rcout << endl;
//     
//     if(iterP < nNewTimes){ // deal with prediction times after the last event
//       
//       if(newtimes[iterP]>etimeMax[iterI]){ // was the computation complete for this event
//         pred_CIF(iterI,iterP) = NA_REAL;  
//       }
//       if(iterP < nNewTimes-1){
//         for(int iterPP = iterP+1; iterPP<nNewTimes ; iterPP++){
// 	  if(newtimes[iterPP]>etimeMax[iterI]){
//             pred_CIF(iterI,iterPP) = NA_REAL;  
//           }else{
//             pred_CIF(iterI,iterPP) = pred_CIF(iterI,iterPP-1);  
//           }
//         }
//       }
//       
//     }
//     
//     if(R_IsNA(t0) == false){ // before t0 fill with NA
//       iterP = 0;
//       while(iterP < nNewTimes && newtimes[iterP]<t0){
//         pred_CIF(iterI,iterP) = NA_REAL;
//         iterP++;
//       }
//      
//     }
//     
//   }
//   
//   return(pred_CIF);
// }
