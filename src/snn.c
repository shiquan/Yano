#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>
#include <assert.h>

SEXP knn2snn(SEXP _knn, SEXP prune)
{
    switch(TYPEOF(_knn)) {
    case INTSXP:
        break;
    default:
        error("'knn' must be an integer matrix");
    }

    PROTECT(_knn);

    int ncell = nrows(_knn);
    int k = ncols(_knn);
    double prune_score = asReal(prune);

    int *knn = INTEGER(_knn);

    int n = 0, m = 1000;
    int *xi = (int*) R_Calloc(ncell + 1, int);
    int *xj = malloc(sizeof(int)*m);
    double *xx = malloc(sizeof(double)*m);
    if (xj == NULL || xx == NULL) {
        free(xj); free(xx); R_Free(xi);
        error("Failed to allocate memory for SNN construction.");
    }

    int *labels = (int*) R_Calloc(ncell, int);
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
                    int *new_xj = realloc(xj, m*sizeof(int));
                    double *new_xx = realloc(xx, m*sizeof(double));
                    if (new_xj == NULL || new_xx == NULL) {
                        free(new_xj); free(new_xx);
                        free(xj); free(xx); R_Free(xi); R_Free(labels);
                        error("Failed to reallocate memory for SNN construction.");
                    }
                    xj = new_xj;
                    xx = new_xx;
                }
                xj[n] = j;
                xx[n] = (double)cnt/(k*2.0 -cnt);
                if (xx[n] >= prune_score) 
                    n++;
            }
        }
    }
    xi[i] = n;

    cholmod_common c;
    M_R_cholmod_start(&c);
        
    CHM_SP ans = M_cholmod_allocate_sparse(ncell, ncell, n, TRUE, TRUE, 0, CHOLMOD_REAL, &c);

    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;

    memcpy(ap, xi, sizeof(int)*(ncell+1));
    memcpy(ai, xj, sizeof(int)*n);
    memcpy(ax, xx, sizeof(double)*n);

    free(xj);
    free(xx);
    R_Free(xi);
    R_Free(labels);
    M_cholmod_finish(&c);

    UNPROTECT(1);
    
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "N", R_NilValue);    
}
