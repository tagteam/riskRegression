#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>
#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <R.h>
#include <Rinternals.h>

#define ME(matrix,row,col) (((matrix)->entries)[(col) * ((matrix)->nr) + (row)])
#define ME3(matrix3,dim1,row,col) (((matrix3)->entries)[(dim1)*(((matrix3)->nr)*((matrix3)->nc))+(col)*((matrix3)->nr)+(row)])
#define VE(vector,index) (((vector)->entries)[(index)])
#define oops(s) {error((s));}
#define max(a,b) ( ((a) > (b)) ? (a) : (b) )
#define min(a,b) ( ((a) > (b)) ? (b) : (a) )
#define malloc_mat(NR, NC, M) { (M) = Calloc(1,matrix); ((M)->nr) = (NR); ((M)->nc) = (NC); ((M)->entries) = Calloc(((NR)*(NC)) , double);}
#define malloc_mat3(DIM,NR, NC, M3) {(M3) = Calloc(1,matrix); ((M3)->dim)=(DIM); ((M3)->nr) = (NR); ((M3)->nc) = (NC); ((M3)->entries) = Calloc(((DIM)*(NR)*(NC)) , double);}
#define malloc_vec(L, V) { (V) = Calloc(1,vector); ((V)->length) = (L); ((V)->entries) = Calloc((L), double);}

typedef struct{
  int dim;
  int nr;
  int nc;
  double *entries;
} matrix3;

typedef struct{
  int nr;
  int nc;
  double *entries;
} matrix;

typedef struct{
  int length;
  double *entries;
} vector;

typedef struct{
  double timec;
  int callc;
} counter;

/* void malloc_mat(int *nrow, int *ncol, matrix *M); */

void free_mat3(matrix3 *M);

void free_mat(matrix *M);

/* void malloc_vec(int *length, vector *V); */

void free_vec(vector *V);

int nrow_matrix(matrix *M);

int ncol_matrix(matrix *M);

int length_vector(vector *v);

void print_a_matrix(matrix *M);

extern void F77_SUB(dpotri)(const char* uplo, const int* n,
			    double* a, const int* lda, int* info);

extern void F77_SUB(dpotrf)(const char* uplo, const int* n,
			    double* a, const int* lda, int* info);

extern void F77_SUB(dgemm)(const char *transa, const char *transb, const int *m, \
			   const int *n, const int *k, const double *alpha,\
			   const double *a, const int *lda,\
			   const double *b, const int *ldb,\
			   const double *beta, double *c, const int *ldc);

extern void F77_SUB(dgemv)(const char *trans, const int *m, const int *n,
			   const double *alpha, const double *a, const int *lda,
			   const double *x, const int *incx, const double *beta,
			   double *y, const int *incy);

extern void F77_SUB(dgetrf)(const int* m, const int* n, double* a, const int* lda,
			    int* ipiv, int* info);

extern void F77_SUB(dgetri)(const int* n, double* a, const int* lda,
                 int* ipiv, double* work, const int* lwork, int* info);

extern void F77_SUB(dqrdc2)(double *x, int *ldx, int *n, int *p,
			    double *tol, int *rank,
			    double *qraux, int *pivot, double *work);

extern void F77_SUB(dtrco)(double*, int*, int*, double*, double*, int*);

extern void F77_SUB(dgecon)(const char* norm, const int* n,
			    const double* a, const int* lda,
			    const double* anorm, double* rcond,
			    double* work, int* iwork, int* info);

extern double F77_NAME(dlange)(const char* norm, const int* m, const int* n,
			       const double* a, const int* lda, double* work);


void MtM(matrix *M, matrix *A);


void cumsumM(matrix *M, matrix *Mout,int rev,int weighted,double *weights); 
void cumsumM1pM2(matrix *M1, matrix *M2,matrix *At[],int rev,int weighted,double *weights,int nindex, int *index);
void cumsumMpM(matrix *M1,matrix *At[],int rev,int weighted,double *weights, int nindex,int *index);

void invertSPD(matrix *A, matrix *AI);
void invertSPDunsafe(matrix *A, matrix *AI);

void Mv(matrix *M, vector *v1, vector *v2);

void vM(matrix *M, vector *v1, vector *v2);

vector *vec_star(vector *v1, vector *v2, vector *v3);
  
double vec_sum(vector *v);

double vec_prod(vector *v1,vector *v2);

double vec_min(vector *v, int *imin);
  
void mat_zeros(matrix *M);
  
void vec_zeros(vector *v);
 
void print_mat(matrix *M);

void print_vec(vector *v);
 
vector *extract_row(matrix *M, int row_to_get, vector *v);

void replace_row(matrix *M, int row_to_set, vector *v);

void vec_add(vector *v1, vector *v2, vector *v3);

vector *scl_vec_mult(double scalar, vector *v1, vector *v2);

matrix *scl_mat_mult(double scalar, matrix *m1, matrix *m2);

matrix *mat_copy(matrix *m1, matrix *m2);

vector *vec_copy(vector *v1, vector *v2);

void mat_subsec(matrix *m1, int rowStart, int colStart,
		       int rowStop, int colStop, matrix *m2);

matrix *mat_transp(matrix *m1, matrix *m2);

void vec_subtr(vector *v1, vector *v2, vector *v3);

void mat_subtr(matrix *m1, matrix *m2, matrix *m3);

void mat_add(matrix *m1, matrix *m2, matrix *m3);

void vec_add_mult(vector *v1, vector *v2, double s, vector *v3);

void MtA(matrix *M, matrix *A, matrix *Mout);

void MAt(matrix *M, matrix *A, matrix *Mout);

void invert(matrix *A, matrix *AI);
void invertS(matrix *A, matrix *AI,int silent);
void invertUnsafe(matrix *A, matrix *AI);
void invertUnsafeS(matrix *A, matrix *AI,int silent);

void cholesky(matrix *A, matrix *AI);
void choleskyunsafe(matrix *A, matrix *AI);

void MxA(matrix *M, matrix *A, matrix *Mout);

void R_CheckUserInterrupt(void);

void print_clock(clock_t *intime,int i);

void update_clock(clock_t *intime, counter *C);

void zcntr(counter *C);

void print_counter(int i, counter *C);

void head_matrix(matrix *M);

void head_vector(vector *V);

void identity_matrix(matrix *M);

void malloc_mats(int nrow, int ncol, ...);

void malloc_vecs(int length, ...);

void free_mats(matrix **M, ...);

void free_vecs(vector **V, ...);

vector *vec_ones(vector *v);

void replace_col(matrix *M, int col_to_set, vector *v);

vector *extract_col(matrix *M, int col_to_get, vector *v);

void Cpred(double *cum,int *nx,int *px,double *xval,int *nxval,double *pred); 

void sindex(int *index, double *jump, double *eval, int *njump, int *neval,int *strict);

void nclusters(int *npers,int *clusters, int *nclust, int *mclust);

void clusterindex(int *clusters, int *nclust,int *npers,int *idclust, int *clustsize,
	          int *mednum,int *num,int *firstclustid); 


void clusterindexdata(int *clusters,int *nclust,int *npers,int *idclust,int *clustsize,int *mednum,
		int *num,double *data, int *p,double *nydata); 
 
void comptest(double *times,int *Ntimes,int *px,double *cu,double *vcu,
	      double *vcudif,int *antsim,double *test,double *testOBS,
	      double *Ut,double *simUt,matrix **W4t,int *weighted,int *antpers); 

double tukey(double x,double b);

void smoothB(double *designX,int *nx,int *p,double *bhat,int *nb,double *b,int *degree,int *coef);

void comptestfunc(double *times,int *Ntimes,int *px,double *cu,double *vcu,
		  double *vcudif,int *antsim,double *test,double *testOBS,
		  double *Ut,double *simUt,matrix **W4t,int *weighted,int *antpers,
		  double *gamma,int *line,double *timepowtest); 

void itfitsemi(double *times,int *Ntimes,double *x,int *delta,int *cause,double *KMc,
	       double *z,int *antpers,int *px,int *Nit,double *score,double *hess,
	       double *est,double *var,int *sim,int *antsim,double *test,
	       double *testOBS,double *Ut,double *simUt,int *weighted, double *gamma,
	       double *vargamma,int *semi,double *zsem,int *pg,int *trans,double *gamma2,
	       int *CA,int *line,int *detail,double *biid,double *gamiid,int *resample,
	       double *timepow,int *clusters,int *antclust,double *timepowtest,int *silent,double *convc,
	       double *weight,double *entry,double *trunkp,int *estimator,int *fixgamma ,int *stratum,
               int *ordertime,int *conservative);
/* , double *ssf); */


void LevenbergMarquardt(matrix *S,matrix *SI,vector *U,vector *delta,double *lm,double *step); 
 
void readXt2(int *antpers,int *nx,int *p,double *designX,
             double *start,double *stop,int *status,int pers,matrix *X,double time);

void readXt(int *antpers,int *nx,int *p,double *designX,double *start,double *stop,
	int *status,int pers,matrix *X,double time,int *clusters,int *cluster,int *id) ;

void readXZt(int *antpers,int *nx,int *px,double *designX,int *pg,double *designG,
		double *start,double *stop,int *status,int pers,matrix *X,
		matrix *WX,matrix *Z,matrix *WZ,double time,int *clusters,
		int *cluster,int *ls,int stat,int l,int *id,int s,int medw);

void readXZtsimple(int *antpers,int *nx,int *px,double *designX,int *pg,double *designG,
		double *start,double *stop,int *status,int pers,matrix *X,
	       	matrix *Z,double time, int s, int *id); 
