// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// C++ functions for calculating the asymptotic covariance matrix for the delongtest function.

// calculates midrank
// [[Rcpp::export]]
vec calculateMidrank(vec z){
  int m = z.size();
  vec wtemp = sort(z);
  uvec index = sort_index(z);
  vec w(z.size()+1);
  for (int i = 0; i <m;i++ ){
    w[i]=wtemp[i];
  }
  w[m]=wtemp[m-1]+1;
  vec t(m,fill::zeros);
  int i = 0;
  int a, b, j;
  while (i < m){
    a = i;
    j = a;
      while (w[j]==w[a]){
        j+=1;
      }
      b = j - 1;
    for (int k = a; k <= b; k++){
      t[k] = (double)(a+b)/2;
    }
    i = b+1;
  }
  int k;
  vec tk(m,fill::zeros);
  for (i = 0; i < m; i++){
    k = index[i];
    tk[k] = t[i]+1;
  }
  return tk;
}

// Fast implementation of the calculation of the covariance matrix.
// Number of rows is the number of observations, for X and Y respectively
// Number of columns is the number of experiments
// note that many of the copying methods may be slow and will need to be optimized
// we note here that pointers and references should be used in order
// to make efficient use of the memory
// [[Rcpp::export]]

NumericMatrix calculateDelongCovarianceFast(NumericMatrix Xs, NumericMatrix Ys){
  int m = Xs.nrow();
  int n = Ys.nrow();

  assert(Xs.ncol()==Ys.ncol());
  int k = Xs.ncol();
  mat V10(k,m);
  mat V01(k,n);
  vec theta(k,fill::zeros);
  for (int r = 0; r < k; r++){
    //theta[r] = 0.0;
    // Make them into armadillo vectors; might be an inefficient and superfluous operation
    // strangely enough we cannot write Xr = as<vec>(Xs(_,r))
    NumericVector Xx = Xs(_,r);
    NumericVector Yy = Ys(_,r);
    vec Xr = as<vec>(Xx);
    vec Yr = as<vec>(Yy);

    //Rcout << "Xr" << Xr << "\n\n\n";
    /*std::cout << "Yr\n";
    */
    //Rcout << "Yr" << Yr << "\n\n\n";
    //Yr.print();
    //std::cout << "\n\n\n\n\n\n\n";*/
    // concatenate
    vec Zr = join_cols(Xr,Yr);
    //Rcout << "zr" << Zr << "\n\n\n";
    // calculate midranks
    vec TZr = calculateMidrank(Zr);
    //Rcout << "TZr" << TZr << "\n\n\n";
    //std::cout << "tzr\n";
    //TZr.print();
    vec TXr = calculateMidrank(Xr);
    //Rcout << "TXr" << TXr << "\n\n\n";
    vec TYr = calculateMidrank(Yr);
    //Rcout << "TYr" << TXr << "\n\n\n";
    //std::cout << "tyr\n";
    //TYr.print();
    //std::cout << "\n\n\n\n\n\n\n";
    for (int i = 0; i < m; i++){
      // FIgure out how Xr is ordered compared to Zr ???
      //?1=TZr[i];
      //
      V10(r,i)=(TZr[i]-TXr[i])/((double) n);
      theta[r]+=TZr[i];
    }
    //Rcout << "r: " << r << "\n" << "V10" << V10.t() << "\n\n\n";
    theta[r]=theta[r]/(double (m*n))-(double) (m+1)/(2*n);
    //Rcout << "r, theta" << r <<  theta;
    for (int j = 0; j < n; j++){
      // FIgure out how Yr is ordered compared to Zr ???
      //?2=TZr[m+j];
      V01(r,j)=1.0-(TZr[j+m]-TYr[j])/((double) m);
    }
    //Rcout << "r: " << r << "\n" <<"V01" << V01.t()  << "\n\n\n";
  }
  mat S(k,k);
  // can be done more effectively
  /*for (int r = 0; r < k; r++){
    for (int s = 0; s < k; s++){
      double s10, s01 = 0.0;
      for (int i = 0; i < m ; i++){
        s10 += (V10(r,i)-theta[r])*(V10(s,i)-theta[s]);
      }
      Rcout << "s10 " << s10 << "\n";
      for (int j = 0; j < n; j++){
        s01 += (V01(r,j)-theta[r])*(V01(s,j)-theta[s]);
      }
      Rcout << "s01 " << s01 << "\n";
      S(r,s)=s10/((double) (m-1))+ s01/(double(n-1));
      Rcout << "result" << S(r,s) << "\n";
      // inprincple also S(s,r) could be specified here
    }
  }
  Rcout << "S is: " << S;*/

  mat s10 = arma::cov(V10.t());
  //Rcout << "s10" << s10 << "\n";
  mat s01 = arma::cov(V01.t());
  //Rcout << "s01" << s01 << "\n";
  //Rcout << "m is " << m << "\n";
  //Rcout << "n is " << n << "\n";
  S = s01/((double) n)+s10/((double) m);
  return wrap(S);
}


// Function is useful in the case where dolist has length 0 and we don't want the entire covariance matrix
// [[Rcpp::export]]
NumericVector calculateDelongDiagonal(int nauc, int nCases, int nControls, NumericMatrix tmn, NumericMatrix tmp){
  //std::cout << nauc << "\n\n";;
  NumericVector SEOfS(nauc);
  for (int j = 0; j < nauc; j++){
    vec temp1(nCases,fill::zeros);
    for (int r = 0; r < nCases; r++){
      for (int k = 0; k < nControls; k++){
        if (tmn(j,k) < tmp(j,r)) {
          temp1[r]+=1.0;
        }
        else if (tmn(j,k) == tmp(j,r)) {
          temp1[r]+=0.5;
        }
      }
    }
    //std::cout << "Temp1: " << temp1.n_rows << "\n";
    //temp1.print();
    //std::cout << "\n\n\n\n\n\n\n";

    vec temp2(nControls,fill::zeros);
    for (int r = 0;  r < nControls; r++){
      for (int k = 0; k < nCases; k++){
        if (tmp(j,k) > tmn(j,r)) {
          temp2[r]+=1.0;
        }
        else if (tmp(j,k) == tmn(j,r)) {
          temp2[r]+=0.5;
        }
      }
    }
    //std::cout << arma::var(temp1);

    //temp2.print();
    /*std::cout << "Temp2: ";
    temp1.print();
    std::cout << "\n";
    std::cout << "length of it: " << temp1.n_rows;*/

    SEOfS[j] = sqrt(arma::var(temp1)/(nCases*nControls*nControls) + arma::var(temp2)/(nCases*nCases*nControls));
    //std::cout << "SE: " << SEOfS[j] << "\n";
  }
  return wrap(SEOfS);
}

/*Pure rewrite of
for (i in 1:nCases) {
    V10[i, ] <- rowSums(tmn < tmp[, i]) + 0.5 * rowSums(tmn == tmp[, i])
}*/
// [[Rcpp::export]]
NumericMatrix rowSumsAlt1(NumericMatrix V, NumericMatrix tmn, NumericMatrix tmp) {
  for (int r = 0; r < V.nrow(); r++) {
    uvec ans(V.ncol(),fill::zeros);
    for (int j = 0; j < tmn.ncol(); j++) {
      for (int i = 0; i < tmn.nrow(); i++) {
        if (tmn(i,j)<tmp(i,r)) {
          ans[i] += 1.0;
        }
        else if (tmn(i,j) == tmp(i,r)) {
          ans[i] += 0.5;
        }
      }
    }
    for (int j = 0; j < V.ncol();j++){
      V(r,j) = ans(j);
    }
  }
  return wrap(V);
}

/*Pure rewrite of
 for (i in 1:nControls) {
    V01[i, ] <- rowSums(tmp > tmn[, i]) + 0.5 * rowSums(tmp == tmn[, i])
 }*/
// [[Rcpp::export]]
NumericMatrix rowSumsAlt2(NumericMatrix V, NumericMatrix tmn, NumericMatrix tmp) {
  for (int r = 0; r < V.nrow(); r++) {
    uvec ans(V.ncol(),fill::zeros);
    for (int j = 0; j < tmp.ncol(); j++) {
      for (int i = 0; i < tmp.nrow(); i++) {
        if (tmp(i,j)>tmn(i,r)) {
          ans[i] += 1.0;
        }
        else if (tmp(i,j) == tmn(i,r)) {
          ans[i] += 0.5;
        }
      }
    }
    for (int j = 0; j < V.ncol();j++){
      V(r,j) = ans(j);
    }
  }
  return wrap(V);
}

