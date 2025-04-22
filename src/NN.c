#ifdef _OPENMP
#include<omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include "Matrix.h"
#include "cholmod.h"
#include <assert.h>

extern double R_euclidean(double *x, int nr, int nc, int i1, int i2);

static cholmod_common c;

struct val {
    int idx;
    double val;
};
int cmpfunc(const void * a, const void * b)
{
    return ((*(struct val*)b).val - (*(struct val*)a).val);
}
int cmpfunc2(const void * a, const void * b)
{
    return ((*(struct val*)b).idx - (*(struct val*)a).idx);
}

// distance between columns in sparseMatrix
SEXP query_nn(SEXP ref, SEXP query, SEXP k_nn, SEXP _threads, SEXP _ndim)
{
    int n_thread = asInteger(_threads);
    int knn = asInteger(k_nn);
    int ndim = asInteger(_ndim);
    int nref = nrows(ref);
    int ncol = ncols(ref);
    
    int nquery = length(query);

    if (knn > nquery) {
        knn = nquery;
    }
    if (ncol < ndim) {
        ndim = ncol;
    }

    R_xlen_t nz;
    nz = (R_xlen_t)nref*knn;
    if (nz > INT_MAX) return(mkString("Too large size."));

    CHM_SP ans = M_cholmod_allocate_sparse(nref, nref, nz, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) return(mkString("Failed to create sparse matrix"));
    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;
    int i;
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < nref; ++i) {
        int n = i * knn;
        int j;
        struct val qval[nquery];
        for (j = 0; j < nquery; ++j) {
            int idx = INTEGER(query)[j]-1;
            qval[j].idx = j;
            qval[j].val = R_euclidean(REAL(ref), nref, ncol, i, idx);
        }
        qsort(qval, nquery, sizeof(struct val), cmpfunc);
        ap[i] = n;
        for (j = knn; j < nquery; ++j) {
            qval[j].idx = INT_MAX;
        }
        qsort(qval, nquery, sizeof(struct val), cmpfunc2);
        for (j = 0; j < knn; ++j) {
            ai[n] = qval[j].idx;
            ax[n] = qval[j].val;
            n++;
        }

        // ap[i+1] = n;        
    }
    
    ap[nref] = nz;
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "N", R_NilValue);
}


