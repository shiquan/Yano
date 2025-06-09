#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>

static cholmod_common c;

SEXP imputation1(SEXP _x, SEXP idx, SEXP _W)
{
    CHM_SP x = AS_CHM_SP__(_x);
    CHM_SP W = AS_CHM_SP__(_W);
    if (W->stype) return mkString("W cannot be symmetric");

    CHM_SP t = M_cholmod_transpose(x, (int)x->xtype, &c);

    const int * tp = (int*)t->p;
    const int * ti = (int*)t->i;
    const double *tx = (double*)t->x;

    const int * xp = (int*)x->p;
    const int * xi = (int*)x->i;
    const double * xx = (double*)x->x;

    const int * wp = (int*)W->p;
    const int * wi = (int*)W->i;
    const double * wx = (double*)W->x;

    R_CheckUserInterrupt();

    int ncells = W->nrow;
    int nfeature = x->nrow;
    
    int n = 0, m = 10000;
    int *yi = malloc(m*sizeof(int));
    double *yx = malloc(m*sizeof(double));
    int yp[length(idx)+1];
    
    int i;
    double w0[ncells];
    int nl = length(idx);
    for (i = 0; i < nl; ++i) { // foreach cell
        yp[i] = n;
        int ii = INTEGER(idx)[i];
        memset(w0, 0, sizeof(double)*ncells);
        int j, k;
        // init weight array
        for (k = wp[ii-1]; k < wp[ii]; ++k) {
            w0[wi[k]] = wx[k];
        }
        
        // foreach feature
        for (j = 0; j < nfeature; ++j) {
            double new = 0;
            for (k = tp[j]; k < tp[j+1]; ++k) {
                int idx = ti[k];
                if (w0[idx] > 0) {
                    new+=w0[idx] * tx[k];                    
                }
            }
            if (new == 0) continue;
            if (n == m) {
                m = m*2;
                yi = realloc(yi, sizeof(int)*m);
                yx = realloc(yx, sizeof(double)*m);
            }
            yi[n] = j;
            yx[n] = new;
            n++;
        }
    }
    yp[i] = n;
    
    M_cholmod_free_sparse(&t, &c);
    CHM_SP ans = M_cholmod_allocate_sparse(nfeature, nl, n, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) return(mkString("Failed to create sparse matrix"));
    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;
    memcpy(ap, yp, sizeof(int)*(nl+1));
    memcpy(ai, yi, sizeof(int)*n);
    memcpy(ax, yx, sizeof(double)*n);
    free(yi);
    free(yx);

    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "N", R_NilValue);
}
