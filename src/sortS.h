// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

int compareS (const void * a, const void * b);

struct structS {
  double time;
  int status;
  double eXb;
  int index;
};

////// sort functions

void sortS(vector<double>& time, vector<int>& status, vector<double>& eXb, uvec& index, int n){
//bool sortS(NumericVector& time, IntegerVector& status, NumericVector& eXb, int n){
  // warning time, status and eXb are called by reference and modified during the execution of the function
  
  // definitions
  structS * dataS = NULL;
 
  // initialisation
  dataS = (structS*) malloc(n * sizeof(structS));
  // if (dataS==NULL) {
  //   return(false);
  // }
  
  for (int i = 0; i < n; ++i) {
    dataS[i].time = time[i];
    dataS[i].status = status[i];
    dataS[i].eXb = eXb[i];
    dataS[i].index = index[i];
  }
  
  // qsort
  qsort(dataS, n, sizeof(structS), compareS);
  
  // export 
  for (int i = 0; i < n; ++i) {
    time[i] = dataS[i].time;
    status[i] = dataS[i].status;
    eXb[i] = dataS[i].eXb;
    index[i] =  dataS[i].index;
  }
  
  free(dataS);
  
  // return(true);
}

////// to feed qsort
int compareS (const void * a, const void * b){ // sort first by time and in case of equality by status (censored last)
  if ((*(structS*)a).time<(*(structS*)b).time) return -1;
  if ((*(structS*)a).time==(*(structS*)b).time){
    if ((*(structS*)a).status == 0){
      if ((*(structS*)b).status == 0) return 0; return -1;
    }else{
      if ((*(structS*)b).status == 0) return -1; return 0;
    }
  } 
  if ((*(structS*)a).time>(*(structS*)b).time) return 1;
  return NA_REAL;
}
