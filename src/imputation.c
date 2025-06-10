#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>

static cholmod_common c;

int scatter(CHM_SP A, int j, double beta, int *w, double *x, int mark, int *yi, int nz)
{
    int i, p;
    int *Ap = A->p;
    int *Ai = A->i;
    double *Ax = A->x;
    for (p = Ap[j]; p < Ap[j+1]; ++p) {
        i = Ai[p];
        if (w[i] < mark) { // init w and x
            w[i] = mark;
            yi[nz++] = i;
            x[i] = beta*Ax[p];
        } else {
            x[i] += beta*Ax[p];
        }
    }
    return nz;
}
SEXP imputation1(SEXP _x, SEXP idx, SEXP _W)
{
    CHM_SP x = AS_CHM_SP__(_x);
    CHM_SP W = AS_CHM_SP__(_W);
    if (W->stype) return mkString("W cannot be symmetric");

    const int * xp = (int*)x->p;
    const int * xi = (int*)x->i;
    const double * xx = (double*)x->x;

    const int * wp = (int*)W->p;
    const int * wi = (int*)W->i;
    const double * wx = (double*)W->x;

    R_CheckUserInterrupt();

    int ncells = W->nrow;
    int nfeature = x->nrow;
    
    uint64_t n = 0, m = 10000;
    if (m < nfeature) m = nfeature;
    int *yi = malloc(m*sizeof(int));
    double *yx = malloc(m*sizeof(double));
    
    int nl = length(idx);

    int yp[nl+1];    
    int w0[nfeature];
    double x0[nfeature];
    memset(w0, 0, sizeof(int)*nfeature);
    memset(x0, 0, sizeof(double)*nfeature);

    int i;
    for (i = 0; i < nl; ++i) { // foreach cell
        yp[i] = n;
        int ii = INTEGER(idx)[i]-1;
        int j;
        for (j = wp[ii]; j < wp[ii+1]; ++j) {
            if (n >= m - nfeature) {
                m = m*2;
                if (m < 0) error("Failed to allocate data.");
                yi = realloc(yi, sizeof(int)*m);
                yx = realloc(yx, sizeof(double)*m);

                if (yi == NULL || yx == NULL) error("Failed to allocate data.");
            }
            n = scatter(x, wi[j], wx[j], w0, x0, ii+1, yi, n);
        }
        for (j = yp[i]; j < n; ++j) yx[j] = x0[yi[j]];
    }
    yp[i] = n;
    
    CHM_SP ans = M_cholmod_allocate_sparse(nfeature, nl, n, FALSE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) return(mkString("Failed to create sparse matrix"));
    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;
    memcpy(ap, yp, sizeof(int)*(nl+1));
    memcpy(ai, yi, sizeof(int)*n);
    memcpy(ax, yx, sizeof(double)*n);
    free(yi);
    free(yx);

    M_cholmod_sort(ans, &c);
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "N", R_NilValue);
}
