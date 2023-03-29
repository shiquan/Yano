#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

void reg(const char* trans, const int* m, const int* n,
         const int* nrhs, double* a, const int* lda,
         double* b, const int* ldb,
         double* work, const int* lwork, int* info)
{
    F77_CALL(dgels)(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info);
}

SEXP mat_mult_blas(SEXP _A, SEXP _a_rows, SEXP _B, SEXP _b_cols, SEXP _k)
{
    char * TRANSA = "N";
    char * TRANSB = "N";
    double * A = REAL(_A);
    double * B = REAL(_B);
    int K = asInteger(_k);
    int M = asInteger(_a_rows);
    int N = asInteger(_b_cols);
    int LDA = K;
    int LDB = N;
    double ALPHA = 1.0;
    double BETA = 0.0;
    SEXP _result;
    PROTECT(_result = allocMatrix(REALSXP, M, N));
    double * result = REAL(_result);  F77_CALL(dgemm)(
        TRANSA, TRANSB, &M, &N, &K, &ALPHA,
        A, &LDA, B, &LDB, &BETA, result, &N);  UNPROTECT(1);
    return(_result);
}
