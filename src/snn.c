#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>
#include <assert.h>

static cholmod_common c;

SEXP knn2snn(SEXP _knn, SEXP prune)
{
    switch(TYPEOF(_knn)) {
    case LGLSXP:
    case INTSXP:
    case REALSXP:
        break;
    default:
        error("'knn' must be numeric");
    }

    if(TYPEOF(_knn) != INTSXP) _knn = coerceVector(_knn, INTSXP);
    PROTECT(_knn);

    int ncell = nrows(_knn);
    int k = ncols(_knn);
    double prune_score = asReal(prune);

    int *knn = INTEGER(_knn);

    int n = 0, m = 1000;
    int xi[ncell+1];
    int *xj = malloc(sizeof(int)*m);
    double *xx = malloc(sizeof(double)*m);
    
    int labels[ncell];
    int i, j;
    int l;
    for (i = 0; i < ncell; ++i) {
        memset(labels, 0, sizeof(int)*ncell);
        xi[i] = n;
        for (l = 0; l < k; ++l) {
            int idx = knn[ncell*l+i]-1;
            labels[idx] = 1;
        }
        for (j = i +1; j < ncell; ++j) {
            int cnt = 0;
            for (l = 0; l < k; ++l) {
                int idx = knn[ncell*l+j]-1;
                if (labels[idx]) cnt++;
            }
            if (cnt) {
                if (n == m) {
                    m = m *2;
                    xj = realloc(xj, m*sizeof(int));
                    xx = realloc(xx, m*sizeof(double));
                }
                xj[n] = j;
                xx[n] = (double)cnt/(k*2.0 -cnt);
                if (xx[n] >= prune_score) 
                    n++;
            }
        }
    }
    xi[i] = n;

    CHM_SP ans = M_cholmod_allocate_sparse(ncell, ncell, n, TRUE, TRUE, 0, CHOLMOD_REAL, &c);

    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;

    memcpy(ap, xi, sizeof(int)*(ncell+1));
    memcpy(ai, xj, sizeof(int)*n);
    memcpy(ax, xx, sizeof(double)*n);

    free(xj);
    free(xx);

    UNPROTECT(1);
    
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "N", R_NilValue);    
}
