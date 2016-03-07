#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

int compareS (const void * a, const void * b);
int compareSS (const void * a, const void * b);

struct structS {
  double time;
  int status;
  double Xb;
};

////// sort functions

bool sortS(vector<double>& time, vector<int>& status, vector<double>& Xb, int n){
//bool sortS(NumericVector& time, IntegerVector& status, NumericVector& Xb, int n){
  // warning time, status and Xb are called by reference and modified during the execution of the function
  
  // definitions
  structS * dataS = NULL;
 
  // initialisation
  dataS = (structS*) malloc(n * sizeof(structS));
  if (dataS==NULL) {
    return(false);
  }
  
  for (int i = 0; i < n; ++i) {
    dataS[i].time = time[i];
    dataS[i].status = status[i];
    dataS[i].Xb = Xb[i];
  }
  
  // qsort
  qsort(dataS, n, sizeof(structS), compareS);
  
  // export 
  for (int i = 0; i < n; ++i) {
    time[i] = dataS[i].time;
    status[i] = dataS[i].status;
    Xb[i] = dataS[i].Xb;
  }
  
  free(dataS);
  
  return(true);
}


////// to feed qsort
int compareS (const void * a, const void * b){
  if ((*(structS*)a).time<(*(structS*)b).time) return -1;
  if ((*(structS*)a).time==(*(structS*)b).time) return 0;
  if ((*(structS*)a).time>(*(structS*)b).time) return 1;
  return NA_REAL;
}

