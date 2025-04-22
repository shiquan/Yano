#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>

#include "utils.h"

// # nocov end
static cholmod_common c;

extern double R_euclidean(double *x, int nr, int nc, int i1, int i2);

SEXP refinedDist(SEXP _dist, SEXP coord, SEXP _dd)
{
    double dd = asReal(_dd);

    switch(TYPEOF(coord)) {
    case LGLSXP:
    case INTSXP:
    case REALSXP:
        break;
    default:
        error("'coord' must be numeric");
    }

    if(TYPEOF(coord) != REALSXP) coord = coerceVector(coord, REALSXP);
    PROTECT(coord);
    int nr = nrows(coord);
    int nc = ncols(coord);

    CHM_SP dist = AS_CHM_SP__(_dist);
    
    int *dp = (int*)dist->p;
    int *di = (int*)dist->i;
    double *dx = (double*)dist->x;

    int ncell = dist->ncol;
    CHM_SP ans = M_cholmod_copy_sparse(dist, &c);
    double *ax = (double*)ans->x;
    int i;
    int j;
    for (i = 0; i < ncell; ++i) {
        if (dp[i] == dp[i+1]) {
            error("No neighbors for cell %d", i+1);
        }
        int idx;
        for (j = dp[i]; j < dp[i+1]; ++j) {
            int i0 = di[j];
            if (i == i0) { //self identity
                dx[j] = 0.0;
            } else {
                double d = R_euclidean(REAL(coord), nr, nc, i, i0);
                d = d/dd;
                //ax[j] = dx[j]*d*d;
                ax[j] = d;
            }
        }
    }
    UNPROTECT(1);
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "N", R_NilValue);
}

SEXP minIdx(SEXP x, SEXP inv)
{
    CHM_SP m = AS_CHM_SP__(x);
    Rboolean invert = asLogical(inv);

    int *ap = (int*)m->p;
    int *ai = (int*)m->i;
    double *ax = (double*)m->x;

    int nc = m->ncol;
    SEXP ival = PROTECT(allocVector(INTSXP, nc));
    
    int i;
    int j;
    for (i = 0; i < nc; ++i) {
        if (ap[i] == ap[i+1]) {
            error("No neighbors for cell %d", i+1);
        }
        double d = -1;
        int idx;
        for (j = ap[i]; j < ap[i+1]; ++j) {
            if (d < 0) {
                d = ax[j];
                idx = ai[j];
            }
            if (d == 0.0) break; // continue; // for self identity
            
            if (invert) {
                if (d < ax[j]) {
                    d = ax[j];
                    idx = ai[j];
                }
            } else {
                if (d > ax[j]) {
                    d = ax[j];
                    idx = ai[j];
                }
            }
        }
        INTEGER(ival)[i] = idx+1;
    }

    UNPROTECT(1);
    return ival;
}
