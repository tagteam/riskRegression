#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void itfit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

extern SEXP riskRegression_baseHaz_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP riskRegression_calcE_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP riskRegression_colCenter_cpp(SEXP, SEXP);
extern SEXP riskRegression_colCumSum(SEXP);
extern SEXP riskRegression_colMultiply_cpp(SEXP, SEXP);
extern SEXP riskRegression_colScale_cpp(SEXP, SEXP);
extern SEXP riskRegression_colSumsCrossprod(SEXP, SEXP, SEXP);
extern SEXP riskRegression_calcSeHazard_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP riskRegression_calcSeCif_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP riskRegression_IFbeta_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);                                     
extern SEXP riskRegression_IFbetaApprox_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP riskRegression_IFlambda0_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP riskRegression_predictCIF_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP riskRegression_quantileProcess_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP riskRegression_rowCenter_cpp(SEXP, SEXP);
extern SEXP riskRegression_rowCumSum(SEXP);
extern SEXP riskRegression_rowMultiply_cpp(SEXP, SEXP);
extern SEXP riskRegression_rowScale_cpp(SEXP, SEXP);
extern SEXP riskRegression_rowSumsCrossprod(SEXP, SEXP, SEXP);
extern SEXP riskRegression_sliceMultiply_cpp(SEXP, SEXP);
extern SEXP riskRegression_sliceScale_cpp(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"itfit", (DL_FUNC) &itfit, 50},
    {NULL, NULL, 0}
};


static const R_CallMethodDef CallEntries[] = {
    {"riskRegression_baseHaz_cpp",         (DL_FUNC) &riskRegression_baseHaz_cpp,         10},
    {"riskRegression_calcE_cpp",           (DL_FUNC) &riskRegression_calcE_cpp,            6},
    {"riskRegression_colCenter_cpp",       (DL_FUNC) &riskRegression_colCenter_cpp,        2},
    {"riskRegression_colCumSum",           (DL_FUNC) &riskRegression_colCumSum,            1},
    {"riskRegression_colMultiply_cpp",     (DL_FUNC) &riskRegression_colMultiply_cpp,      2},
    {"riskRegression_colScale_cpp",        (DL_FUNC) &riskRegression_colScale_cpp,         2},
    {"riskRegression_colSumsCrossprod",    (DL_FUNC) &riskRegression_colSumsCrossprod,     3},
    {"riskRegression_calcSeHazard_cpp",    (DL_FUNC) &riskRegression_calcSeHazard_cpp,    26},
    {"riskRegression_calcSeCif_cpp",       (DL_FUNC) &riskRegression_calcSeCif_cpp,       33},
    {"riskRegression_IFbeta_cpp",          (DL_FUNC) &riskRegression_IFbeta_cpp,          10},
    {"riskRegression_IFbetaApprox_cpp",    (DL_FUNC) &riskRegression_IFbetaApprox_cpp,     6},
    {"riskRegression_IFlambda0_cpp",       (DL_FUNC) &riskRegression_IFlambda0_cpp,       16},   
    {"riskRegression_predictCIF_cpp",      (DL_FUNC) &riskRegression_predictCIF_cpp,      15},
    {"riskRegression_quantileProcess_cpp", (DL_FUNC) &riskRegression_quantileProcess_cpp,  6},
    {"riskRegression_rowCenter_cpp",       (DL_FUNC) &riskRegression_rowCenter_cpp,        2},
    {"riskRegression_rowCumSum",           (DL_FUNC) &riskRegression_rowCumSum,            1},
    {"riskRegression_rowMultiply_cpp",     (DL_FUNC) &riskRegression_rowMultiply_cpp,      2},
    {"riskRegression_rowScale_cpp",        (DL_FUNC) &riskRegression_rowScale_cpp,         2},
    {"riskRegression_rowSumsCrossprod",    (DL_FUNC) &riskRegression_rowSumsCrossprod,     3},
    {"riskRegression_sliceMultiply_cpp",   (DL_FUNC) &riskRegression_sliceMultiply_cpp,    2},
    {"riskRegression_sliceScale_cpp",      (DL_FUNC) &riskRegression_sliceScale_cpp,       2},
    {NULL, NULL, 0}
};

void R_init_riskRegression(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
