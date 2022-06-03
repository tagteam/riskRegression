// [[Rcpp::depends(RcppArmadillo)]]
#include "arma-wrap.h"
using namespace Rcpp;

// Same as colMeans(A*b) for a matrrix A and vector b
// [[Rcpp::export]]
NumericVector columnMeanWeight(NumericMatrix A, NumericVector x){
  int nrows = A.nrow();
  int ncols = A.ncol();
  NumericVector ans(ncols);
  for (int j = 0; j < ncols; j++){
    double sum = 0.0;
    for (int i = 0; i < nrows; i++){
      sum += A(i,j)*x[i];
    }
    ans[j] = sum/nrows;
  }
  return ans;
}


// Same as colSums(b*(A+1)) for a matrix A and vector b
// [[Rcpp::export]]
NumericVector T3CalculationHelper(NumericVector x, NumericMatrix A){
  int nrows = A.nrow();
  int ncols = A.ncol();
  NumericVector ans(ncols);
  for (int j = 0; j < ncols; j++){
    double sum = 0.0;
    for (int i = 0; i < nrows; i++){
      sum += (A(i,j)+1.0)*x[i];
    }
    ans[j] = sum;
  }
  return ans;
}

// [[Rcpp::export]]
NumericMatrix htijCalculationHelper(NumericVector mcase, NumericVector mcontrol,NumericVector wcase,NumericVector wcontrol,int n, int nrows, int ncols){
  NumericMatrix ans(nrows,ncols);
  for (int j = 0; j < ncols; j++){
    for (int i = 0; i < nrows; i++){
      if (mcase(i) > mcontrol(j)) {
        ans(i,j) = wcase(i)*wcontrol(j)*n*n;
      }
      else if (mcase(i) == mcontrol(j)) {
        ans(i,j) = 0.5* wcase(i)*wcontrol(j)*n*n;
      }
      else {
        ans(i,j) = 0.0;
      }
    }
  }
  return ans;
}

// [[Rcpp::export]]
NumericMatrix rowSumsCrossprodSpec(arma::mat &X, arma::mat &Y){
  return(wrap(arma::sum(X,1).t()*(Y+1)));
}

// [[Rcpp::export]]
NumericMatrix colSumsCrossprodSpec(arma::mat &X, arma::mat &Y){
  return(wrap(arma::sum(X,0)*(Y+1)));
}

// arma::mat htijCalculation(arma::vec &mcase, arma::vec &mcontrol,arma::vec &wcase,arma::vec &wcontrol,int n, int nrows, int ncols){
//   arma::mat ans(nrows,ncols);
//   for (int j = 0; j < ncols; j++){
//     for (int i = 0; i < nrows; i++){
//       if (mcase(i) > mcontrol(j)) {
//         ans(i,j) = wcase(i)*wcontrol(j)*n*n;
//       }
//       else if (mcase(i) == mcontrol(j)) {
//         ans(i,j) = 0.5* wcase(i)*wcontrol(j)*n*n;
//       }
//     }
//   }
//   return ans;
// }

// arma::vec Cpp_colSums(const arma::mat& x) {
//   int nr = x.n_rows, nc = x.n_cols;
//   arma::vec ans(nc);
//   for (int j = 0; j < nc; j++) {
//     double sum = 0.0;
//     for (int i = 0; i < nr; i++) {
//       sum += x(i, j);
//     }
//     ans[j] = sum;
//   }
//   return ans;
// }

// arma::vec Cpp_rowSums(const arma::mat& x) {
//   int nr = x.n_rows, nc = x.n_cols;
//   arma::vec ans(nr);
//   for (int j = 0; j < nc; j++) {
//     for (int i = 0; i < nr; i++) {
//       ans[i] += x(i, j);
//     }
//   }
//   return ans;
// }

// arma::vec sindexHelp(arma::vec &jumpTimes, arma::vec &evalTimes ) {
//   int neval = evalTimes.n_elem;
//   int n = jumpTimes.n_elem;
//   arma::vec sortedJump = sort(jumpTimes);
//   arma::vec sortedEval = sort(evalTimes);
//   arma::vec ans(neval);
//   int i,t;
//   ans[0] = 0;
//   i = 0;
//   for (t=0;t<neval;t++){
//     while(i<n && sortedJump[i]<=sortedEval[t]) i++;
//     ans[t] = i;
//   }
//   int k;
//   arma::vec ans2(neval);
//   arma::uvec finInd = sort_index(sort_index(evalTimes));
//   for (int i = 0; i < neval; i++){
//     k = finInd[i];
//     ans2[i] = ans[k];
//   }
//   return ans2;
// }

// arma::mat T3Calculation(double &hathtstar, double &F01t, arma::vec &vectTisupt, arma::vec &fi1t, arma::mat &MC){
//   int nrows = MC.n_rows;
//   int ncols = MC.n_cols;
//   arma::vec ans(ncols);
//   for (int j = 0; j < ncols; j++){
//     double sum = 0.0;
//     for (int i = 0; i < nrows; i++){
//       sum += (MC(i,j)+1.0)*fi1t[i];
//     }
//     ans[j] = sum;
//   }
//   //return ans;
//   return hathtstar*(sum(vectTisupt) + 1/F01t * ans-nrows);
// }

// // [[Rcpp::export]]
// NumericVector getInfluenceCurveAUCSurvival(int n, arma::vec &time, arma::vec &risk, arma::vec &cases, arma::vec &controls, arma::vec &ipcwControls, arma::vec &ipcwCases, arma::mat &MC){
//   int nbCases = arma::accu(cases);
//   int nbControls = n-nbCases;
//   arma::vec uniqueRisk = unique(risk);
//   if (nbCases==0 || nbControls ==0 || uniqueRisk.n_elem==1) {
//     return(wrap(arma::vec(n)));
//   } 
//   double F01t = arma::accu(ipcwCases);
//   double St = arma::accu(ipcwControls);
//   arma::vec mcase(nbCases);
//   arma::vec wcase(nbCases);
//   arma::vec mcontrols(nbControls);
//   arma::vec wcontrols(nbControls);
//   arma::vec fi1t(n);
//   arma::vec timeCases(nbCases);
//   int caseCounter = 0;
//   int controlCounter = 0;
//   for (int i = 0; i < n; i++){
//     if (cases[i]) {
//       mcase[caseCounter] = risk[i];
//       wcase[caseCounter] = ipcwCases[i];
//       timeCases[caseCounter] = time[i];
//       fi1t[i] = ipcwCases[i]*n;
//       caseCounter++;
//     }
//     else {
//       mcontrols[controlCounter] = risk[i];
//       wcontrols[controlCounter] = ipcwControls[i];
//       controlCounter++;
//     }
//   }
//   // Rcout << mcase << "\n";
//   // Rcout << wcase << "\n";
//   // Rcout << mcontrols << "\n";
//   // Rcout << wcontrols << "\n";
//   arma::mat htij1 = htijCalculation(mcase,mcontrols,wcase,wcontrols,n,nbCases,nbControls);
//   //return wrap(htij1(0,0));
//   arma::vec colSumshtij1(n);
//   arma::vec rowSumshtij1(n);
//   arma::vec* temp1 = new arma::vec(nbCases);
//   *temp1 = Cpp_rowSums(htij1);
//   arma::vec* temp2 = new arma::vec(nbControls);
//   *temp2 = Cpp_colSums(htij1);
//   caseCounter = 0;
//   controlCounter = 0;
//   for (int i = 0; i < n; i++){
//     if (cases[i]) {
//       colSumshtij1[i] = (*temp1)[caseCounter];
//       caseCounter++;
//     }
//     else {
//       rowSumshtij1[i] = (*temp2)[controlCounter];
//       controlCounter++;
//     }
//   }
//   delete temp1;
//   delete temp2;
//   arma::vec vectTisupt = double(n)/double(nbControls)*controls;
//   double hathtstar = accu(htij1)/(n*n);
//   arma::vec uniqueTime = unique(time);
//   arma::vec sind = sindexHelp(uniqueTime,timeCases);
//   arma::mat MCTiCases(timeCases.n_elem,MC.n_cols);
//   caseCounter = 0;
//   int ind;
//   for (int i = 0; i < timeCases.n_elem; i++){
//     ind = sind[i]-1;
//     for (int j = 0; j < MC.n_cols; j++){
//       MCTiCases(i,j) = MC(ind,j)+1.0;
//     } 
//   }  
//   arma::rowvec T1 = 1.0/n*arma::sum(htij1,1).t()*MCTiCases;
//   arma::vec T3 = T3Calculation(hathtstar,F01t,vectTisupt,fi1t,MC);
//   // forget ties for now
//   arma::vec termijak(T3.n_elem);
//   for (int i = 0; i < T3.n_elem;i++){
//     termijak[i] = (T1[i]-T3[i])/(F01t*St);
//   }
//   arma::vec termikaj = (rowSumshtij1 - n*hathtstar)/(F01t*St);
//   arma::vec termjkai = (colSumshtij1 - n*hathtstar*(vectTisupt+(1/F01t)*(fi1t-F01t)))/(F01t*St);
//   return wrap((termijak+termikaj+termjkai)/n);
// }

// // [[Rcpp::export]]
// int sindexSingle(arma::vec &jumpTimes, double &evalTime ) {
//   arma::vec sortedJump = sort(jumpTimes);
//   int i = 0;
//   while(i<jumpTimes.n_elem && sortedJump[i]<=evalTime) i++;
//   return(i);
// }
// 
// arma::vec columnMeanWeight(arma::mat &A, arma::vec x){
//   int nrows = A.n_rows;
//   int ncols = A.n_cols;
//   arma::vec ans(ncols);
//   for (int j = 0; j < ncols; j++){
//     double sum = 0.0;
//     for (int i = 0; i < nrows; i++){
//       sum += A(i,j)*x[i];
//     }
//     ans[j] = sum/nrows;
//   }
//   return ans;
// }
// 
// // [[Rcpp::export]]
// NumericMatrix getInfluenceCurveBrier(double &t, arma::vec &time,arma::vec &ic0,arma::vec  &residuals,arma::vec &Wti, arma::vec &Wt,arma::mat &ICG, String censModel, int nthTimes){
//   int n = residuals.n_elem;
//   if (censModel == "cox"){
//     return 0;
//   }
//   else {
//     if (ICG.n_rows==0){
//       return wrap(residuals-mean(residuals)); 
//     }
//     else {
//       double Brier = mean(residuals);
//       arma::vec hit1(residuals.n_elem);
//       arma::vec hit2(residuals.n_elem);
//       for (int i = 0; i< residuals.n_elem; i++){
//         if (time[i] <= t){
//           hit2[i] = residuals[i];
//         }
//         else {
//           hit1[i] = residuals[i];
//         }
//       }
//       arma::vec uniqueTime = unique(time);
//       int ind = sindexSingle(uniqueTime,t);
//       if (ind > 0){
//         arma::vec Int0tdMCsurEffARisk(ICG.n_cols);
//         for (int i = 0; i < ICG.n_cols; i++){
//           Int0tdMCsurEffARisk[i] = ICG(ind-1,i);
//         }
//         return(wrap(residuals-Brier+mean(hit1)*Int0tdMCsurEffARisk+columnMeanWeight(ICG,hit2)));
//       }
//       else {
//         return(wrap(residuals-Brier+columnMeanWeight(ICG,hit2)));
//       }
//     }
//   }
// }
// 
