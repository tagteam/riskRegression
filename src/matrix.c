#include <stdlib.h>
//#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <R_ext/Random.h>
#include "matrix.h"

void free_mat(matrix *M){
  Free(M->entries);
  Free(M);
}

void free_mat3(matrix3 *M){
  Free(M->entries);
  Free(M);
}

void free_vec(vector *V){
  Free(V->entries);
  Free(V);
}

int nrow_matrix(matrix *M){
  return M->nr;
}

int ncol_matrix(matrix *M){

  return M->nc;

}

int length_vector(vector *v){

  return v->length;

}

void print_a_matrix(matrix *M){

  int j, k;
  for(j=0; j < nrow_matrix(M); j++){
    for(k = 0; k < ncol_matrix(M); k++){
      Rprintf("%+7.7g ", ME(M,j,k));
    }
    Rprintf("\n");
  }  

}

/* DPOTRI - compute the inverse of a real symmetric positive */
/* definite matrix A using the Cholesky factorization A = U**T*U */
/* or A = L*L**T computed by DPOTRF */
extern void F77_SUB(dpotri)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);


/* DPOTRF - compute the Cholesky factorization of a real */
/* symmetric positive definite matrix A */
extern void F77_SUB(dpotrf)(const char* uplo, const int* n,
		 double* a, const int* lda, int* info);


extern void F77_SUB(dgemm)(const char *transa, const char *transb, const int *m,
		const int *n, const int *k, const double *alpha,
		const double *a, const int *lda,
		const double *b, const int *ldb,
		const double *beta, double *c, const int *ldc);

/* DGEMV - perform one of the matrix-vector operations */
/* y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y,  */

extern void F77_SUB(dgemv)(const char *trans, const int *m, const int *n,
		const double *alpha, const double *a, const int *lda,
		const double *x, const int *incx, const double *beta,
		double *y, const int *incy);


/* DGETRF - compute an LU factorization of a general M-by-N */
/* matrix A using partial pivoting with row interchanges */
extern void
F77_SUB(dgetrf)(const int* m, const int* n, double* a, const int* lda,
                 int* ipiv, int* info);


/* DGETRI - compute the inverse of a matrix using the LU */
/* factorization computed by DGETRF */
extern void
F77_SUB(dgetri)(const int* n, double* a, const int* lda,
                 int* ipiv, double* work, const int* lwork, int* info);


// cumsum of matrix apply(X,2,cusum)
// rev=1 apply(X[n:1,],2,cumsum)[n:1,]
// for rev=1 possible to return only apply(X[n:1,],2,cumsum)[nindex,]
void cumsumM(matrix *M, matrix *Mout,int rev,int weighted,double *weights)   {
  int i,j,p=ncol_matrix(M),n=nrow_matrix(M); 
  double lweights[n]; 

  matrix *temp; 
  malloc_mat(n,p,temp); 

  if( !( 
	(ncol_matrix(M) == ncol_matrix(Mout)) ))  {
    oops("Error: dimensions in cumsumM\n");
  }

  for(i=0;i<n;i++)if (weighted==0)lweights[i]=1.0; else lweights[i]=weights[i]; 

  if (rev==0) {
    for(j = 0; j < n; j++) ME(Mout,0,j)=ME(M,0,j)*lweights[0]; 
    for(i = 1; i < n; i++) 
      for(j = 0; j < n; j++) ME(Mout,i,j)= ME(Mout,i-1,j)+ ME(M,i,j)*lweights[i]; 
  }

  if (rev==1) {
    matrix *temp; 
    malloc_mat(n,p,temp); 

    for(j = 0; j < p; j++) ME(temp,0,j)=ME(M,n-1,j)*lweights[n-1]; 
    for(i = 1; i < n; i++) 
      for(j = 0; j < p; j++) ME(temp,i,j)= ME(temp,i-1,j)+ME(M,n-1-i,j)*lweights[n-1-i]; 

    for(i = 0; i < n; i++) 
      for(j = 0; j < p; j++) ME(Mout,i,j)=ME(temp,n-1-i,j);
    free_mat(temp); 
  }

}

// returns cumsum X^T %*% Z  , rev=1 does the reverse sum 
// and index only returns certain indeces of the cumsum only for rev=1
// see a cumsumM
void cumsumM1pM2(matrix *M1, matrix *M2,matrix *At[],int rev,int weighted,double *weights,int nindex, int *index){
  int i,j,k,p1=ncol_matrix(M1),p2=ncol_matrix(M2),n=nrow_matrix(M1); 
  double lweights[n]; 

  if( !( (p1== nrow_matrix(At[0])) && 
	 (p2== ncol_matrix(At[0])) )){
    oops("Error: dimensions in cumsumM1pM2\n");
  }

  for(i = 0; i < n; i++)if (weighted==0)lweights[i]=1.0; else lweights[i]=weights[i]; 

  if (rev==0) {
    for(j = 0; j < p1; j++) 
      for(k = 0; k < p2; k++) ME(At[0],j,k)=ME(M1,0,j)*ME(M2,0,k)*lweights[0]; 

    for(i = 1; i < n; i++) 
      for(j = 0; j < p1; j++) 
	for(k = 0; k < p2; k++) ME(At[i],j,k)= ME(At[i-1],j,k)+
	  ME(M1,i,j)*ME(M2,i,k)*lweights[i]; 
  }

  if (rev==1) {
  matrix *temp[n],*temp1[n]; 
  for(i = 0; i < n; i++) { 
	  malloc_mat(p1,p2,temp[i]); malloc_mat(p1,p2,temp1[i]); }

    for(j = 0; j < p1; j++) 
      for(k = 0; k < p2; k++) ME(temp[0],j,k)=ME(M1,n-1,j)*ME(M2,n-1,k)*lweights[n-1]; 

    for(i = 1; i < n; i++) 
      for(j = 0; j < p1; j++) 
	for(k = 0; k < p2; k++) ME(temp[i],j,k)= ME(temp[i-1],j,k)+
	  ME(M1,n-i-1,j)*ME(M2,n-i-1,k)*lweights[n-i-1]; 

  for(i = 0; i < n; i++)  mat_copy(temp[i],temp1[n-i-1]);  

if (nindex>0) 
  for(i = 0;i<nindex;i++)  {
	  mat_copy(temp1[index[i]],At[i]);  
  }
  else for(i = 0; i < n; i++)  mat_copy(temp[i],At[n-i-1]);  

  for(i = 0;i<n;i++)  { free_mat(temp[i]); free_mat(temp1[i]); }  
}

}

// returns cumsum X^T %*% Z  , rev=1 does the reverse sum 
// and index only returns certain indeces of the cumsum only for rev=1
// see a cumsumM
void cumsumMpM(matrix *M1,matrix *At[],int rev,int weighted,double *weights,
		int nindex,int *index){
int i,j,k,p1=ncol_matrix(M1),p2=ncol_matrix(M1),n=nrow_matrix(M1); 
double lweights[n]; 

// head_matrix(M1); head_matrix(At[0]); 

  if( !( 
	(p1== nrow_matrix(At[0])) && 
	(p2== ncol_matrix(At[0])) )){
    oops("Error: dimensions in cumsumMpM\n");
  }

  for(i = 0; i < n; i++)if (weighted==0)lweights[i]=1.0; else lweights[i]=weights[i]; 

if (rev==0) {
  for(j = 0; j < p1; j++) 
  for(k = 0; k < p2; k++) ME(At[0],j,k)=ME(M1,0,j)*ME(M1,0,k)*lweights[0]; 

  for(i = 1; i < n; i++) 
  for(j = 0; j < p1; j++) 
  for(k = 0; k < p2; k++) ME(At[i],j,k)= ME(At[i-1],j,k)+
                          ME(M1,i,j)*ME(M1,i,k)*lweights[i]; 
}

if (rev==1) {
  matrix *temp[n],*temp1[n]; 
  for(i = 0; i < n; i++) { 
	  malloc_mat(p1,p2,temp[i]); malloc_mat(p1,p2,temp1[i]); }

  for(j = 0; j < p1; j++) 
  for(k = 0; k < p2; k++) ME(temp[0],j,k)=ME(M1,n-1,j)*ME(M1,n-1,k)*lweights[n-1]; 

  for(i = 1; i < n; i++) 
  for(j = 0; j < p1; j++) 
  for(k = 0; k < p2; k++) ME(temp[i],j,k)= ME(temp[i-1],j,k)+
	                  ME(M1,n-i-1,j)*ME(M1,n-i-1,k)*lweights[n-i-1]; 

  for(i = 0; i < n; i++)  mat_copy(temp[i],temp1[n-i-1]);  

if (nindex>0) 
  for(i = 0;i<nindex;i++)  {
	  mat_copy(temp1[index[i]],At[i]);  
  }
  else for(i = 0; i < n; i++)  mat_copy(temp[i],At[n-i-1]);  

  for(i = 0;i<n;i++)  { free_mat(temp[i]); free_mat(temp1[i]); }  
}

}


// Performs A := t(M) %*% M, where A is an nRowM x nColM matrix, 
// and A is an nColM x nColM matrix
void MtM(matrix *M, matrix *A){

  char transa = 't';
  char transb = 'n';
  double alpha = 1.0;
  double beta = 0.0;
  int m = ncol_matrix(M);
  int n = ncol_matrix(M);
  int k = nrow_matrix(M);
  int lda = nrow_matrix(M);
  int ldb = nrow_matrix(M);
  int ldc = ncol_matrix(M);

  if( !(nrow_matrix(A) == ncol_matrix(M) && 
	ncol_matrix(A) == ncol_matrix(M)) ){
    oops("Error: dimensions in MtM\n");
  }

  // Ensure that M and A do not occupy the same memory. 
  if(M != A){
    // the results of 1.0 * t(M) %*% M + 0.0 * c is stored in c
    F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &alpha, M->entries, &lda,
		    M->entries, &ldb, &beta, A->entries, &ldc);
  } else {
    // if M and A occupy the same memory, store the results in a
    // temporary matrix. 
    matrix *temp;
    malloc_mat(nrow_matrix(A),ncol_matrix(A),temp);

    F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &alpha, M->entries, &lda,
		    M->entries, &ldb, &beta, temp->entries, &ldc);

    // Copy these results into A, then remove the temporary matrix
    mat_copy(temp,A);

    free_mat(temp);
  }
}

// Does cholesky of := A, where A is symmetric positive definite, of order *n
void cholesky(matrix *A, matrix *AI){ // {{{ 

  if( !(nrow_matrix(A)  == ncol_matrix(A) && 
	nrow_matrix(AI) == ncol_matrix(AI) &&
	nrow_matrix(A)  == ncol_matrix(AI)) ){
    oops("Error: dimensions in invertSPD\n");
  }

  // Ensure that A and AI do not occupy the same memory. 
  if(A != AI){
//	  printf(" er her\n"); 
    choleskyunsafe(A, AI);
  } else {
    // if M and A occupy the same memory, store the results in a
    // temporary matrix. 
    matrix *temp;
    malloc_mat(nrow_matrix(AI),ncol_matrix(AI),temp);

    choleskyunsafe(A, temp);
    
    // Copy these results into AI, then remove the temporary matrix
    mat_copy(temp,AI);    
    free_mat(temp);
  }

} // }}} 

// cholesky := A, where A is symmetric positive definite, of order *n
void choleskyunsafe(matrix *A, matrix *AI){ // {{{ 
  //unsafe because it assumes A and AI are both square and of the same
  //dimensions, and that they occupy different memory

//  char uplo = 'U'; // lower version 
  int i, j;
  int n = nrow_matrix(A);
//  int lda = n; // matrix A has dimensions *n x *n
  int info = -999;
//    double rcond;
//  int pivot[n];
//  double z[n];
//  double qraux[n];
//  double work[2*n];
//  int rank = 0;
//  int job=1;
//  double tol = 1.0e-07;
  
// First copy the matrix A into the matrix AI
//  print_mat(A); 
   mat_copy(A,AI); 
//  print_mat(AI); 

//  printf("sssssssssss======================\n"); 
//  job = 1; // Indicates that AI is upper triangular
//  rcond = 999.0;

    // First find the Cholesky factorization of A,
    // stored as an upper triangular matrix
    char uplo1 = 'U'; // lower version 
    F77_CALL(dpotrf)(&uplo1, &n, AI->entries, &n, &info); 

    // Lastly turn the vector a into the matrix AI
    // Take only the lower triangular portion, since this 
    // is the relevant part returned by dpotrf
    for(i = 0; i < n; i++){
      for(j = 0; j < i; j++){
	      ME(AI,i,j) = 0; 
      }
    }
//    print_mat(AI); 

//    Rprintf("in chol \n"); 
//    printf("======================\n"); 
//    print_mat(A); 
//    print_mat(AI); 
//    printf(" check chol back\n"); 
//    matrix *tmp; 
//    malloc_mat(n,n,tmp); 
//    MtM(AI,tmp); 
//    print_mat(tmp); 
//    free_mat(tmp); 

} // }}} 

// Does AI := inverse(A), where A is symmetric positive definite, of order *n
void invertSPD(matrix *A, matrix *AI){

  if( !(nrow_matrix(A)  == ncol_matrix(A) && 
	nrow_matrix(AI) == ncol_matrix(AI) &&
	nrow_matrix(A)  == ncol_matrix(AI)) ){
    oops("Error: dimensions in invertSPD\n");
  }

  // Ensure that A and AI do not occupy the same memory. 
  if(A != AI){
    invertSPDunsafe(A, AI);
  } else {
    // if M and A occupy the same memory, store the results in a
    // temporary matrix. 
    matrix *temp;
    malloc_mat(nrow_matrix(AI),ncol_matrix(AI),temp);

    invertSPDunsafe(A, temp);
    
    // Copy these results into AI, then remove the temporary matrix
    mat_copy(temp,AI);    
    free_mat(temp);
  }

}

// Does AI := inverse(A), where A is symmetric positive definite, of order *n
void invertSPDunsafe(matrix *A, matrix *AI){
  //unsafe because it assumes A and AI are both square and of the same
  //dimensions, and that they occupy different memory

  char uplo = 'U';
  int i, j;
  int n = nrow_matrix(A);
  int lda = n; // matrix A has dimensions *n x *n
  int info = -999;
  double rcond;
  int pivot[n];
  double z[n];
  double qraux[n];
  double work[2*n];
  int rank = 0;
  int job=1;
  double tol = 1.0e-07;
  
  // First copy the matrix A into the matrix AI
  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      ME(AI,i,j) = ME(A,i,j);
    }
  }
    
//  dqrdc(x,ldx,n,p, qraux,jpvt,work,job)
//  F77_CALL(dqrdc)(AI->entries, &n, &n, &n, &rank, qraux, pivot, work,job);
//  dqrdc2(x,ldx,n,p,tol,k,qraux,jpvt,work)
  F77_CALL(dqrdc2)(AI->entries, &n, &n, &n, &tol, &rank, qraux, pivot, work);

  for(i = 0; i < n; i++){
    for(j = 0; j < i; j++){
      ME(AI,j,i) = 0.0;
    }
  }

  job = 1; // Indicates that AI is upper triangular
  rcond = 999.0;
  F77_CALL(dtrco)(AI->entries, &n, &n, &rcond, z, &job);
    
  if(rcond < tol){
    Rprintf("Error in invertSPD: estimated condition number = %7.7e\n",1/rcond); 
    
    for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
	ME(AI,i,j) = 0.0;
      }
    }
  } else {

    for(i = 0; i < n; i++){
      pivot[i] = i+1;
      for(j = 0; j < n; j++){
	ME(AI,i,j) = ME(A,i,j);
      }
    }

    // First find the Cholesky factorization of A,
    // stored as an upper triangular matrix
    F77_CALL(dpotrf)(&uplo, &n, AI->entries, &lda, &info); 

    if(info < 0){
      Rprintf("Error in invertSPD: arg %d of DPOTRF\n",-info);
    } else if(info > 0){
      Rprintf("Error in invertSPD: matrix does not appear to be SPD\n");
    } 
       
    // then use this factorization to compute the inverse of A
    F77_CALL(dpotri)(&uplo, &n, AI->entries, &lda, &info);

    if(info != 0){
      Rprintf("Error in invertSPD: DPOTRI returned info = %d \n",info);
    }

    // Lastly turn the vector a into the matrix AI
    // Take only the upper triangular portion, since this 
    // is the relevant part returned by dpotrf
    for(i = 0; i < n; i++){
      for(j = 0; j < i; j++){
	ME(AI,i,j) = ME(AI,j,i);
      }
    }
    
  }

}

// v2 := M %*% v1
// where M has dims (nrow x ncol)
// and v1 has dims  (ncol x  1  )
// amd v2 has dims  (nrow x  1  )
void Mv(matrix *M, vector *v1, vector *v2){

  char trans = 'n';
  double alpha = 1.0;
  double beta = 0.0;
  int incx = 1;
  int incy = 1;
  int nrow = nrow_matrix(M);
  int ncol = ncol_matrix(M);

  if( !(length_vector(v1) == ncol && 
	length_vector(v2) == nrow) ){
    oops("Error: dimensions in Mv\n");
  }
  
  // Ensure that v1 and v2 do not occupy the same memory. 
  if(v1 != v2){
    F77_CALL(dgemv)(&trans, &nrow, &ncol, &alpha, M->entries, &nrow,
		    v1->entries, &incx, &beta, v2->entries, &incy);
  } else {
    // if v1 and v2 occupy the same memory, store the results in a
    // temporary vector. 
    vector *temp;
    malloc_vec(length_vector(v2),temp);

    F77_CALL(dgemv)(&trans, &nrow, &ncol, &alpha, M->entries, &nrow,
		    v1->entries, &incx, &beta, temp->entries, &incy);
    
    // Copy these results into A, then remove the temporary matrix
    vec_copy(temp,v2);    
    free_vec(temp);
  }
  
}

// v2 := v1 %*% matrix
// where v1 has dims     (1    x nrow)
// and matrix has dims   (nrow x ncol)
// amd v2 has dims       (1    x ncol)
void vM(matrix *M, vector *v1, vector *v2){

  char trans = 't';
  double alpha = 1.0;
  double beta = 0.0;
  int incx = 1;
  int incy = 1;
  int nrow = nrow_matrix(M);
  int ncol = ncol_matrix(M);

  if( !(length_vector(v1) == nrow && 
	length_vector(v2) == ncol) ){
    oops("Error: dimensions in vM\n");
  }
  
  // Ensure that v1 and v2 do not occupy the same memory. 
  if(v1 != v2){
    F77_CALL(dgemv)(&trans, &nrow, &ncol, &alpha, M->entries, &nrow,
		    v1->entries, &incx, &beta, v2->entries, &incy);
  } else {
    // if v1 and v2 occupy the same memory, store the results in a
    // temporary vector. 
    vector *temp;
    malloc_vec(length_vector(v2),temp);

    F77_CALL(dgemv)(&trans, &nrow, &ncol, &alpha, M->entries, &nrow,
		    v1->entries, &incx, &beta, temp->entries, &incy);
    
    // Copy these results into A, then remove the temporary matrix
    vec_copy(temp,v2);    
    free_vec(temp);
  }
}

// v3 := v1 * v2, where * is the Hadamard (componentwise) product of the 
// two vectors, which is the same as * does in R for vectors of the same length
vector *vec_star(vector *v1, vector *v2, vector *v3){
  
  int i;
  int n = length_vector(v1);

  if( !(length_vector(v2) == n && 
	length_vector(v3) == n) ){
    oops("Error: dimensions in vec_star\n");
  }

  for(i = 0; i < n; i++){
    VE(v3,i) = VE(v1,i)*VE(v2,i);
  }

  return(v3);
  
}

// := v1^T * v2, inner product  
// two vectors
double vec_prod(vector *v1, vector *v2){
  
  double sum = 0.0;
  int i;
  int n = length_vector(v1);

  if( !(length_vector(v2) == n) ){
    oops("Error: dimensions in vec_star\n");
  }

  for(i = 0; i < n; i++){
    sum += VE(v1,i)*VE(v2,i);
  }

  return sum;
  
}

// Sums the entries of a vector of length n
double vec_sum(vector *v){
  
  double sum = 0.0;
  int i;
  int n = length_vector(v);

  for(i = 0; i < n; i++){
    sum += VE(v,i);
  }
  return sum;
  
}

// Sums the entries of a vector of length n
vector *vec_ones(vector *v){
  
  int i;
  int n = length_vector(v);

  for(i = 0; i < n; i++){
    VE(v,i) = 1.0;
  }

  return(v);

}

// Returns the minimum of the entries of a vector of length n
double vec_min(vector *v, int *imin){
  
  double Min = VE(v,0);
  int i;
  int n = length_vector(v);
  *imin = 0;

  for(i = 1; i < n; i++){
    if(VE(v,i) < Min){
      Min = VE(v,i);
      *imin = i;
    }
  }
  return Min;
  
}


// set all entries of an *nrow x *ncol matrix M to zero
void mat_zeros(matrix *M){
  
  int j, k;
  
  for(j=0; j < nrow_matrix(M); j++){
    for(k = 0; k < ncol_matrix(M); k++){
      ME(M,j,k) = 0.0;
    }
  }  
  
}

// set all entries of vector v of length *length to zero
void vec_zeros(vector *v){
  
  int j;
  
  for(j=0; j < length_vector(v); j++){
    VE(v,j) = 0.0;
  }
  
}

// Simple I/O function that prints a matrix
void print_mat(matrix *M){
 
  int j, k;

  Rprintf("Matrix nrow=%d ncol=%d \n",nrow_matrix(M),ncol_matrix(M)); 
  for(j=0; j < nrow_matrix(M); j++){
    for(k = 0; k < ncol_matrix(M); k++){
//      Rprintf("%5.5g ", ME(M,j,k));
      // Rprintf("%+15.15g ", ME(M,j,k));
      Rprintf("%lf ", ME(M,j,k));
    }
    Rprintf("\n");
  }  
    Rprintf("\n");
}

// Simple I/O function that prints the top of a matrix
void head_matrix(matrix *M){
 
  int j, k;

  Rprintf("head:Matrix nrow=%d ncol=%d \n",nrow_matrix(M),ncol_matrix(M)); 
  for(j=0; j < min(nrow_matrix(M),6); j++){
    for(k = 0; k < min(ncol_matrix(M),6); k++){
      //Rprintf("%5.5g ", ME(M,j,k));
      Rprintf("%lf ", ME(M,j,k));
    }
    Rprintf("\n");
  }  
    Rprintf("\n");
}

// Simple I/O function that prints the first few entries of a vector
void head_vector(vector *V){
 
  int j;

  Rprintf("head:Vector lengthn=%d \n",length_vector(V)); 
  for(j=0; j < min(length_vector(V),6); j++){
    Rprintf("%lf ", VE(V,j));
  }  
  Rprintf("\n");

}


// Simple I/O function that prints a vector
void print_vec(vector *v){
 
  int j;

  Rprintf("Vector lengthn=%d \n",length_vector(v)); 
  for(j=0; j < length_vector(v); j++){
    Rprintf("%lf ", VE(v,j));
  }  
  Rprintf("\n\n");
  
}

// sets v := M[row_to_get,]
vector *extract_row(matrix *M, int row_to_get, vector *v){

  int j;

  if(!(length_vector(v) == ncol_matrix(M))){
    oops("Error: dimensions in extract_row\n");
  }

  if(row_to_get >= 0 && row_to_get < nrow_matrix(M)){ 
    for(j = 0; j < length_vector(v); j++){
      VE(v,j) = ME(M,row_to_get,j);
    }
    return(v);
  } else {
    oops("Error: trying to get an invalid row in 'extract_row'\n");
  }

  return(v);
    
}

// sets M[row_to_get,] := v
void replace_row(matrix *M, int row_to_set, vector *v){

  int j;

  if(!(length_vector(v) == ncol_matrix(M))){
    oops("Error: dimensions in replace_row\n");
  }

  if(row_to_set >= 0 && row_to_set < nrow_matrix(M)){
    for(j = 0; j < ncol_matrix(M); j++){
      ME(M,row_to_set,j) = VE(v,j);
    }
  } else {
    oops("Error: trying to get an invalid row in 'replace_row'\n");
  }
  
}

// v3 := v1 + v2, where the three vectors have length
void vec_add(vector *v1, vector *v2, vector *v3){

  int i;
  int n = length_vector(v1);

  if( !(length_vector(v2) == n && 
	length_vector(v3) == n) ){
    oops("Error: dimensions in vec_addition\n");
  }

  for(i=0; i < n; i++){
    VE(v3,i) = VE(v1,i) + VE(v2,i);
  }

}

// v3 := v1 + s * v2, where the three vectors have length,
// and s is a double scalar
void vec_add_mult(vector *v1, vector *v2, double s, vector *v3){

  int i;
  int n = length_vector(v1);

  if( !(length_vector(v2) == n && 
	length_vector(v3) == n) ){
    oops("Error: dimensions in vec_addition\n");
  }

  for(i=0; i < n; i++){
    VE(v3,i) = VE(v1,i) + s*VE(v2,i);
  }

}


// v2 := scalar * v1, where invec and outvec are vectors of
// length *length, and *scalar is a (double) scalar
vector *scl_vec_mult(double scalar, vector *v1, vector *v2){

  int i;
  int n = length_vector(v1);

  if( !(length_vector(v2) == n) ){
    oops("Error: dimensions in scl_vec_mult\n");
  }

  for(i=0; i < n; i++){
    VE(v2,i) = scalar * VE(v1,i);
  }

  return(v2);
  
}

// m2 := scalar * m1
matrix *scl_mat_mult(double scalar, matrix *m1, matrix *m2){

  int i,j;
  int m = nrow_matrix(m1);
  int n = ncol_matrix(m1);
  
  if( !(nrow_matrix(m1) == m && ncol_matrix(m1) == n) ){
    oops("Error: dimensions in scl_vec_mult\n");
  }

  for(i=0; i < m; i++){
    for(j=0; j < n; j++){
      ME(m2,i,j) = ME(m1,i,j) * scalar;
    }
  }
  
  return(m2);
  
}

// m2 := m1
matrix *mat_copy(matrix *m1, matrix *m2){

  int i,j;
  int m = nrow_matrix(m1);
  int n = ncol_matrix(m1);

  if( !(nrow_matrix(m2) == m && ncol_matrix(m2) == n) ){
    oops("Error: dimensions in copy_matrix\n");
  }

  if(m1 == m2){
    oops("copy_matrix was asked to write one matrix into its own memory\nThere may be an error...\n");
  }

  for(i=0; i < m; i++){
    for(j=0; j < n; j++){
      ME(m2,i,j) = ME(m1,i,j);
    }
  }

  return(m2);
  
}

// v2 := v1
vector *vec_copy(vector *v1, vector *v2){

  int i;
  int l = length_vector(v1);

  if( !(length_vector(v2) == l) ){
    oops("Error: dimensions in copy_vector\n");
  }

  if(v1 == v2){
    oops("copy_vector was asked to write one matrix into its own memory\nThere may be an error...\n");
  }

  for(i=0; i < l; i++){
    VE(v2,i) = VE(v1,i);
  }

  return(v2);
  
}


// m2 := m1
void mat_subsec(matrix *m1, int rowStart, int colStart,
		       int rowStop, int colStop, matrix *m2){

  int i,j;
  int m = nrow_matrix(m1);
  int n = ncol_matrix(m1);

  if( !(nrow_matrix(m2) == (rowStop-rowStart) 
	&& ncol_matrix(m2) == (colStop-colStart)) ){
    oops("Error: dimensions in mat_subsec\n");
  } else if(!(rowStart >= 0 && colStart >= 0
	      && rowStop < m && colStop < n)){
    oops("Error: trying to access non-existing rows or cols in mat_subsec\n");
  }

  if(m1 == m2){
    oops("matrix_subsec was asked to write one matrix into its own memory\nThere may be an error...\n");
  }

  for(i=rowStart; i < rowStop; i++){
    for(j=colStart; j < colStop; j++){
      ME(m2,i-rowStart,j-colStart) = ME(m1,i,j);
    }
  }
  
}

// m2 := t(m1)
matrix *mat_transp(matrix *m1, matrix *m2){

  int i,j;
  int m = nrow_matrix(m1);
  int n = ncol_matrix(m1);

  if( !(ncol_matrix(m2) == m && nrow_matrix(m2) == n) ){
    oops("Error: dimensions in mat_transp\n");
  }

  // Ensure that m1 and m2 do not occupy the same memory. 
  if(m1 != m2){

    for(i=0; i < m; i++){
      for(j=0; j < n; j++){
	ME(m2,j,i) = ME(m1,i,j);
      }
    }

  } else {
    // if v1 and v2 occupy the same memory, store the results in a
    // temporary vector. 
    matrix *temp;
    malloc_mat(nrow_matrix(m2),ncol_matrix(m2),temp);

    for(i=0; i < m; i++){
      for(j=0; j < n; j++){
	ME(temp,j,i) = ME(m1,i,j);
      }
    }
    
    // Copy these results into A, then remove the temporary matrix
    mat_copy(temp,m2);    
    free_mat(temp);
  }

  return(m2);
  
}


// v3 := v1 - v2, where the three vectors have length *length
void vec_subtr(vector *v1, vector *v2, vector *v3){

  int i;
  int n = length_vector(v1);

  if( !(length_vector(v2) == n && 
	length_vector(v3) == n) ){
    oops("Error: dimensions in vec_subtraction\n");
  }

  for(i=0; i < n; i++){
    VE(v3,i) = VE(v1,i) - VE(v2,i);
  }

}

// m3 := m1 - m2, where the three matrix have the same dimentions
void mat_subtr(matrix *m1, matrix *m2, matrix *m3){

  int i,j;
  int m = nrow_matrix(m1);
  int n = ncol_matrix(m1);

  if( !(nrow_matrix(m2) == m && ncol_matrix(m2) == n && 
	nrow_matrix(m3) == m && ncol_matrix(m3) == n) ){
    oops("Error: dimensions in mat_subtr\n");
  }

  for(i=0; i < m; i++){
    for(j=0; j < n; j++){
      ME(m3,i,j) = ME(m1,i,j) - ME(m2,i,j);
    }
  }

}

// m3 := m1 + m2, where the three matrix have the same dimentions
void mat_add(matrix *m1, matrix *m2, matrix *m3){

  int i,j;
  int m = nrow_matrix(m1);
  int n = ncol_matrix(m1);

  if( !(nrow_matrix(m2) == m && ncol_matrix(m2) == n &&
	nrow_matrix(m3) == m && ncol_matrix(m3) == n) ){
    oops("Error: dimensions in mat_subtr\n");
  }

  for(i=0; i < m; i++){
    for(j=0; j < n; j++){
      ME(m3,i,j) = ME(m1,i,j) + ME(m2,i,j);
    }
  }

}

// Performs Mout := t(M) %*% A, where M is an nRowM x nColM matrix, 
// and A is an nRowM x nColA matrix, and Mout is a nColM x nColA matrix
void MtA(matrix *M, matrix *A, matrix *Mout){

  char transa = 't';
  char transb = 'n';
  double alpha = 1.0;
  double beta = 0.0;
  int m = ncol_matrix(M); 
  int n = ncol_matrix(A);
  int k = nrow_matrix(M);
  int lda = nrow_matrix(M);
  int ldb = nrow_matrix(M);
  int ldc = ncol_matrix(M);

  if( !(nrow_matrix(M)    == nrow_matrix(A) && 
	nrow_matrix(Mout) == ncol_matrix(M) &&	
	ncol_matrix(Mout) == ncol_matrix(A)) ){
    oops("Error: dimensions in MtA\n");
  }
 
  // Ensure that Mout does not occupy the same memory as M or A 
  if(Mout != A && Mout != M){
    // the results of 1.0 * t(M) %*% A + 0.0 * Mout is stored in Mout
    F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &alpha, M->entries,
		    &lda, A->entries, &ldb, &beta, Mout->entries, &ldc);
  } else {
    // if M and A occupy the same memory, store the results in a
    // temporary matrix. 
    matrix *temp;
    malloc_mat(nrow_matrix(Mout),ncol_matrix(Mout),temp);

    // the results of 1.0 * t(M) %*% A + 0.0 * Mout is stored in temp
    F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &alpha, M->entries,
		    &lda, A->entries, &ldb, &beta, temp->entries, &ldc);
    
    // Copy these results into A, then remove the temporary matrix
    mat_copy(temp,Mout);    
    free_mat(temp);
  }


}

// Performs Mout := M %*% t(A), where M is an nRowM x nColM matrix, 
// and A is an nRowA x nColM matrix, and Mout is a nRowM x nRowA matrix
void MAt(matrix *M, matrix *A, matrix *Mout){

  char transa = 'n';
  char transb = 't';
  double alpha = 1.0;
  double beta = 0.0;
  int m = nrow_matrix(M); 
  int n = nrow_matrix(A);
  int k = ncol_matrix(M);
  int lda = nrow_matrix(M);
  int ldb = nrow_matrix(A);
  int ldc = nrow_matrix(Mout);

  if( !(ncol_matrix(M)    == ncol_matrix(A) && 
	nrow_matrix(Mout) == nrow_matrix(M) &&	
	ncol_matrix(Mout) == nrow_matrix(A)) ){
    oops("Error: dimensions in MAt\n");
  }
 
  // Ensure that Mout does not occupy the same memory as M or A 
  if(Mout != A && Mout != M){
    // the results of 1.0 * t(M) %*% A + 0.0 * Mout is stored in Mout
    F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &alpha, M->entries,
		    &lda, A->entries, &ldb, &beta, Mout->entries, &ldc);
  } else {
    // if M and A occupy the same memory, store the results in a
    // temporary matrix. 
    matrix *temp;
    malloc_mat(nrow_matrix(Mout),ncol_matrix(Mout),temp);

    // the results of 1.0 * t(M) %*% A + 0.0 * Mout is stored in temp
    F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &alpha, M->entries,
		    &lda, A->entries, &ldb, &beta, temp->entries, &ldc);
    
    // Copy these results into Mout, then remove the temporary matrix
    mat_copy(temp,Mout);    
    free_mat(temp);
  }


}

// Does Ainv := inverse(A), where A is a square matrix
void invert(matrix *A, matrix *Ainv){

  if( !(nrow_matrix(A)  == ncol_matrix(A) && 
	nrow_matrix(Ainv) == ncol_matrix(Ainv) &&
	nrow_matrix(A)  == ncol_matrix(Ainv)) ){
    oops("Error: dimensions in invert\n");
  }

  // Ensure that A and Ainv do not occupy the same memory. 
  if(A != Ainv){
    invertUnsafe(A, Ainv);
  } else {
    // if A and Ainv occupy the same memory, store the results in a
    // temporary matrix. 
    matrix *temp;
    malloc_mat(nrow_matrix(Ainv),ncol_matrix(Ainv),temp);

    invertUnsafe(A, temp);
    
    // Copy these results into Ainv, then remove the temporary matrix
    mat_copy(temp,Ainv);    
    free_mat(temp);
  }
}

// Does Ainv := inverse(A), where A is a square matrix
void invertS(matrix *A, matrix *Ainv,int silent){

  if( !(nrow_matrix(A)  == ncol_matrix(A) && 
	nrow_matrix(Ainv) == ncol_matrix(Ainv) &&
	nrow_matrix(A)  == ncol_matrix(Ainv))){
    oops("Error: dimensions in invert\n");
  }

  // Ensure that A and Ainv do not occupy the same memory. 
  if(A != Ainv){
    invertUnsafeS(A, Ainv,silent);
  } else {
    // if A and Ainv occupy the same memory, store the results in a
    // temporary matrix. 
    matrix *temp;
    malloc_mat(nrow_matrix(Ainv),ncol_matrix(Ainv),temp);

    invertUnsafeS(A, temp,silent);
    
    // Copy these results into Ainv, then remove the temporary matrix
    mat_copy(temp,Ainv);    
    free_mat(temp);
  }
}


// Does Ainv := inverse(A), where A is a square matrix
void invertUnsafe(matrix *A, matrix *Ainv){
  //unsafe because it assumes A and Ainv are both square and of the
  //same dimensions, and that they occupy different memory

  //char uplo = 'U';
  int i, j;
  int n = nrow_matrix(A);
  int lda = n; // matrix A has dimensions n x n
  int *ipiv = malloc(n * sizeof(int));
  int lwork = n * n;
  int info = -999;
  double anorm = -999.0;
  double rcond = -999.0;
  double tol = 1.0e-07;
  double *dwork = malloc(4 * n * sizeof(double));
  int *iwork = malloc(n * sizeof(int));
  double *work = malloc(n * n * sizeof(double));

  // First turn the matrix A into the vector a

  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      ME(Ainv,i,j) = ME(A,i,j);
    }
  }

  anorm = F77_NAME(dlange)("O", &n, &n, Ainv->entries, &lda, dwork);

  // First find the LU factorization of A,
  // stored as an upper triangular matrix
  F77_CALL(dgetrf)(&n, &n, Ainv->entries, &lda, ipiv, &info);

  if(info != 0){
    //Avoid printing this error message
    Rprintf("2 Error in invert: DGETRF returned info = %d \n",info);
    mat_zeros(Ainv);
    print_mat(Ainv); 
  } else {
  
    for(i = 0; i < n; i++){
      iwork[i]= ipiv[i];
    }
    F77_CALL(dgecon)("O", &n, Ainv->entries, &lda, &anorm, &rcond, dwork, iwork,  &info);
    
    if(info != 0){
      //Avoid printing this error message
      Rprintf("1 Error in invert: DGETRF returned info = %d \n",info);
      mat_zeros(Ainv);
      return;
    } 
    
    if(rcond < tol){
      Rprintf("Error in invert: estimated reciprocal condition number = %7.7e\n",rcond); 
      mat_zeros(Ainv);
      return;
    }

    // then use this factorization to compute the inverse of A
    F77_CALL(dgetri)(&n, Ainv->entries, &lda, ipiv, work, &lwork, &info);

    if(info != 0){
      Rprintf("Error in invert: DPOTRI returned info = %d \n",info);
      mat_zeros(Ainv);
    }

    if (fabs(ME(Ainv,0,0))>99999999999999)  { // TS 23-10
      print_mat(Ainv);
      Rprintf("Inversion, unstable large elements  \n");
      mat_zeros(Ainv);
    }
  }
  
  free(work);
  free(iwork);
  free(dwork);
  free(ipiv);
    
}

// Does Ainv := inverse(A), where A is a square matrix, possibly silent
void invertUnsafeS(matrix *A, matrix *Ainv,int silent){
  //unsafe because it assumes A and Ainv are both square and of the
  //same dimensions, and that they occupy different memory

  //char uplo = 'U';
  int i, j;
  int n = nrow_matrix(A);
  int lda = n; // matrix A has dimensions n x n
  int *ipiv = malloc(n * sizeof(int));
  int lwork = n * n;
  int info = -999;
  double anorm = -999.0;
  double rcond = -999.0;
  double tol = 1.0e-07;
  double *dwork = malloc(4 * n * sizeof(double));
  int *iwork = malloc(n * sizeof(int));
  double *work = malloc(n * n * sizeof(double));

  // First turn the matrix A into the vector a

  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){
      ME(Ainv,i,j) = ME(A,i,j);
    }
  }

  anorm = F77_NAME(dlange)("O", &n, &n, Ainv->entries, &lda, dwork);

  // First find the LU factorization of A,
  // stored as an upper triangular matrix
  F77_CALL(dgetrf)(&n, &n, Ainv->entries, &lda, ipiv, &info);

  if(info != 0){
    //Avoid printing this error message
    mat_zeros(Ainv);
    if (silent==0) Rprintf("3 Error in invert: DGETRF returned info = %d \n",info);
  } else {
  
    for(i = 0; i < n; i++){
      iwork[i]= ipiv[i];
    }
    F77_CALL(dgecon)("O", &n, Ainv->entries, &lda, &anorm, &rcond, dwork, iwork,  &info);
    
    if(info != 0){
      //Avoid printing this error message
      mat_zeros(Ainv);
      free(work); free(iwork); free(dwork); free(ipiv);
      if (silent==0) Rprintf("4 Error in invert: DGETRF returned info = %d \n",info);
      return;
    } 
    
    if(rcond < tol ){
      mat_zeros(Ainv);
      free(work); free(iwork); free(dwork); free(ipiv);
      if (silent==0) Rprintf("Error in invert: estimated reciprocal condition number = %7.7e\n",rcond); 
      return;
    }

    // then use this factorization to compute the inverse of A
    F77_CALL(dgetri)(&n, Ainv->entries, &lda, ipiv, work, &lwork, &info);

    if(info != 0 ){
      mat_zeros(Ainv);
      if (silent==0) Rprintf("Error in invert: DPOTRI returned info = %d \n",info);
    }

    if (fabs(ME(Ainv,0,0))>99999999999999 )  { // TS 23-10
      mat_zeros(Ainv);
      if (silent==0) Rprintf("Inversion, unstable large elements  \n");
    }
  }
  
  free(work); free(iwork); free(dwork); free(ipiv);
}


// Performs Mout := M %*% A, where M is an nRowM x nColM matrix, 
// and A is an nColM x nColA matrix, and Mout is a nRowM x nColA matrix
void MxA(matrix *M, matrix *A, matrix *Mout){

	char transa = 'n';
	char transb = 'n';
	double alpha = 1.0;
	double beta = 0.0;
	int m = nrow_matrix(M);
	int n = ncol_matrix(A);
	int k = ncol_matrix(M);
	int lda = nrow_matrix(M);
	int ldb = ncol_matrix(M);
	int ldc = nrow_matrix(M);

	if( !(ncol_matrix(M)    == nrow_matrix(A) && 
				nrow_matrix(Mout) == nrow_matrix(M) &&	
				ncol_matrix(Mout) == ncol_matrix(A)) ){
		oops("Error: dimensions in MxA\n");
	}

  // Ensure that Mout does not occupy the same memory as M or A 
  if(Mout != A && Mout != M){
    // the results of 1.0 * M %*% A + 0.0 * c is stored in c
    // therfore we do not need to initialise c
    F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &alpha, M->entries, &lda,
		    A->entries, &ldb, &beta, Mout->entries, &ldc);
  } else {
    // if M and A occupy the same memory, store the results in a
    // temporary matrix. 
    matrix *temp;
    malloc_mat(nrow_matrix(Mout),ncol_matrix(Mout),temp);

    // the results of 1.0 * M %*% A + 0.0 * c is stored in c
    // therfore we do not need to initialise c
    F77_CALL(dgemm)(&transa, &transb, &m, &n, &k, &alpha, M->entries, &lda,
		    A->entries, &ldb, &beta, temp->entries, &ldc);
    
    // Copy these results into Mout, then remove the temporary matrix
    mat_copy(temp,Mout);    
    free_mat(temp);
  }


}

void print_clock(clock_t *intime, int i){

  clock_t outtime = clock();
  
  Rprintf("### point %d, time %7.7e\n", i, difftime(outtime,*intime));
  
  *intime = outtime;

}

void update_clock(clock_t *intime, counter *C){

  clock_t outtime = clock();
  
  C->timec += difftime(outtime,*intime);
  C->callc++;
  
  *intime = outtime;

}

void zcntr(counter *C){
  C->timec = 0.0;
  C->callc = 0;
}

void print_counter(int i, counter *C){

  Rprintf("### counter %d, time %7.7g, calls %d\n", i, C->timec, C->callc);
  
}

void identity_matrix(matrix *M){

  int i, j;

  if(nrow_matrix(M) != ncol_matrix(M)){
    oops("Error in identity_matrix: dimenions do not match\n");
  }

  for(i = 0; i < nrow_matrix(M); i++){
    for(j = 0; j < nrow_matrix(M); j++){
      if(i == j){
	ME(M,i,j) = 1.0;
      } else {
	ME(M,i,j) = 0.0;
      }
    }
  }
}

void malloc_mats(int nrow, int ncol, ...){

  va_list argp;
  va_start(argp, ncol);
  matrix **M;

  while((M = va_arg(argp, matrix **))){
    malloc_mat(nrow,ncol,*M);
  }

  va_end(argp);

}

void malloc_vecs(int length, ...){

  va_list argp;
  va_start(argp, length);
  vector **V;

  while((V = va_arg(argp, vector **))){
    malloc_vec(length,*V);
  }

  va_end(argp);

}

void free_mats(matrix **M1, ...){

  va_list argp;
  va_start(argp, M1);
  matrix **M;

  free_mat(*M1);

  while((M = va_arg(argp, matrix **))){
    free_mat(*M);
  }

  va_end(argp);

}

void free_vecs(vector **V1, ...){

  va_list argp;
  va_start(argp, V1);
  vector **V;

  free_vec(*V1);

  while((V = va_arg(argp, vector **))){
    free_vec(*V);
  }

  va_end(argp);

}

// sets v := M[,col_to_get]
vector *extract_col(matrix *M, int col_to_get, vector *v){

  int j;

  if(!(length_vector(v) == nrow_matrix(M))){
    oops("Error: dimensions in extract_col\n");
  }

  if(col_to_get >= 0 && col_to_get < ncol_matrix(M)){ 
    for(j = 0; j < length_vector(v); j++){
      VE(v,j) = ME(M,j,col_to_get);
    }
    return(v);
  } else {
    oops("Error: trying to get an invalid column in 'extract_col'\n");
  }

  return(v);
    
}

// sets M[,col_to_set] := v
void replace_col(matrix *M, int col_to_set, vector *v){

  int j;

  if(!(length_vector(v) == nrow_matrix(M))){
    oops("Error: dimensions in replace_col\n");
  }

  if(col_to_set >= 0 && col_to_set < ncol_matrix(M)){
    for(j = 0; j < nrow_matrix(M); j++){
      ME(M,j,col_to_set) = VE(v,j);
    }
  } else {
    oops("Error: trying to get an invalid column in 'replace_col'\n");
  }
  
}


void LevenbergMarquardt(matrix *S,matrix *SI,vector *U,vector *delta,double *lm,double *step)
{ // {{{
  int i,nrow;
  double ss=0;
  matrix *S2; 

  if(!(length_vector(U) == nrow_matrix(S))){ oops("Error: LM : S and U not consistent\n"); }
  if(!(length_vector(U) == length_vector(delta))){ 
       oops("Error: LM : delta and U not consistent\n"); }

  nrow=length_vector(delta); 
  malloc_mat(nrow,nrow,S2); 
  for (i=0;i<nrow;i++) ss=ss+VE(U,i)*VE(U,i);

  mat_copy(S,S2); 

  if  (ss > *lm ) {
     MxA(S,S,S2);
     for (i=0;i<nrow;i++) ME(S2,i,i)=ME(S2,i,i)+min(VE(U,i)*VE(U,i),100);
     invertS(S2,SI,1); MxA(SI,S,S2); Mv(S2,U,delta);
  } else {
    invertS(S2,SI,1); Mv(SI,U,delta);
  }
  if  (*step>0.0001) scl_vec_mult(*step,delta,delta); 
  free_mat(S2); 
} // }}}

void readXt2(int *antpers,int *nx,int *p,double *designX,
             double *start,double *stop,int *status,int pers,matrix *X,double time)  
{ // {{{
  int j,c,count;

   for (c=0,count=0;((c<*nx) && (count!=*antpers));c++)
   {
	if ((start[c]<time) && (stop[c]>=time)) {
	  for(j=0;j<*p;j++){ ME(X,count,j) = designX[j*(*nx)+c]; }
	  if (time==stop[c] && status[c]==1) { pers=count; }
	  count=count+1; 
        } 
   }
} // }}}

void readXt(int *antpers,int *nx,int *p,double *designX,double *start,double *stop,int *status,int pers,matrix *X,double time,int *clusters,int *cluster,int *id) 
{ // {{{
  int j,c,count;

for (c=0,count=0;((c<*nx) && (count!=*antpers));c++){
      if ((start[c]<time) && (stop[c]>=time)) {
	for(j=0;j<*p;j++) {
	  ME(X,id[c],j) = designX[j*(*nx)+c]; }
	  cluster[id[c]]=clusters[c]; 
	if (time==stop[c] && status[c]==1) { pers=id[c]; }
	count=count+1; } 
    }
} // }}}

void readXZt(int *antpers,int *nx,int *px,double *designX,int *pg,double *designG,
		double *start,double *stop,int *status,int pers,matrix *X,
		matrix *WX,matrix *Z,matrix *WZ,double time,int *clusters,
		int *cluster,int *ls,int stat,int l,int *id,int s,int medw)  
{ // {{{ 
int j,c,count,pmax;

pmax=max(*pg,*px); 

for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	if ((start[c]<time) && (stop[c]>=time)) {
	  cluster[id[c]]=clusters[c]; 
	  for(j=0;j<pmax;j++) {
	    if (j<*px) { ME(X,id[c],j)=designX[j*(*nx)+c]; }
	    if (medw==1) { if (j<*px) { ME(WX,id[c],j) =designX[j*(*nx)+c]; }}
	    if (j<*pg) { ME(Z,id[c],j)=designG[j*(*nx)+c]; }
	    if (medw==1) {if (j<*pg) { ME(WZ,id[c],j)=designG[j*(*nx)+c]; } }
	  }
	  if (time==stop[c] && status[c]==1) {
	    pers=id[c];stat=1;l=l+1; ls[l]=s;
	  }
	  count=count+1; 
	}
      }

} // }}}

void readXZtsimple(int *antpers,int *nx,int *px,double *designX,int *pg,double *designG,
		double *start,double *stop,int *status,int pers,matrix *X,
	       	matrix *Z,double time, int s,int *id)  
{ // {{{ 
int j,c,count,pmax;

pmax=max(*pg,*px); 

for (c=0,count=0;((c<*nx) && (count!=*antpers));c++) {
	if ((start[c]<time) && (stop[c]>=time)) {
	 // cluster[id[c]]=clusters[c]; 
	  for(j=0;j<pmax;j++) {
	    if (j<*px) { ME(X,id[c],j)=designX[j*(*nx)+c]; }
	    if (j<*pg) { ME(Z,id[c],j)=designG[j*(*nx)+c]; }
	  }
	  if (time==stop[c] && status[c]==1) { pers=id[c]; }
	  count=count+1; 
	}
      }

} // }}}

