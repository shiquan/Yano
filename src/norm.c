#include <R.h>
#include <Rinternals.h>

#include "utils.h"
#include "dict.h"
#include "Matrix/Matrix.h"


SEXP lognorm(SEXP _X, SEXP cs, SEXP _sf)
{
    cholmod_common c;
    M_R_cholmod_start(&c);

    double sf = asReal(_sf);
    CHM_SP X = AS_CHM_SP__(_X);

    const int *xp = (int*)X->p;
    const int *xi = (int*)X->i;
    const double *xx = (double*)X->x;

    int nrow = X->nrow;
    int ncol = X->ncol;
    CHM_SP ans = M_cholmod_allocate_sparse(nrow, ncol, xp[ncol], TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) return(mkString("Failed to create sparse matrix"));
    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;

    int i, p, n = 0;
    for (i = 0; i < ncol; ++i) {
        ap[i] = n;
        double cs0 = REAL(cs)[i];
        for (p = xp[i]; p < xp[i+1]; ++p, ++n) {
            ai[p] = xi[p];
            ax[p] = log1p(xx[p]/cs0*sf);
        }
    }
    ap[i] = n;
    M_cholmod_finish(&c);
    
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "N", R_NilValue);
}

SEXP norm_by_row(SEXP A) {
    cholmod_common c;
    M_R_cholmod_start(&c);

    CHM_SP X = AS_CHM_SP__(A);
    CHM_SP B = M_cholmod_copy(X, 0, TRUE, &c);
        
    const int *xp = (int*)X->p;
    const int *xi = (int*)X->i;
    const double *xx = (double*)X->x;

    double *bx = (double*)B->x;
    int nc = X->ncol;
    for (int i = 0; i < nc; ++i) {
        double tot = 0;
        for (int j = xp[i]; j < xp[i+1]; ++j) {
            tot += xx[j];            
        }
        if (tot == 0) continue;
        for (int j = xp[i]; j < xp[i+1]; ++j) {
            bx[j] = bx[j]/tot;
        }
    }

    M_cholmod_finish(&c);
    return M_chm_sparse_to_SEXP(B, 1, 0, 0, "N", R_NilValue); 
}
