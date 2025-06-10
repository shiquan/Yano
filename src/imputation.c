#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>

static cholmod_common c;

size_t scatter(CHM_SP A, size_t j, double beta, size_t *w, double *x, size_t mark, size_t *yi, size_t nz)
{
    size_t i, p;
    size_t *Ap = A->p;
    size_t *Ai = A->i;
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

    const size_t * xp = (size_t*)x->p;
    const size_t * xi = (size_t*)x->i;
    const double * xx = (double*)x->x;

    const size_t * wp = (size_t*)W->p;
    const size_t * wi = (size_t*)W->i;
    const double * wx = (double*)W->x;

    size_t ncells = W->nrow;
    size_t nfeature = x->nrow;
    
    size_t n = 0, m = 10000;
    if (m < nfeature) m = nfeature;
    size_t *yi = malloc(m*sizeof(size_t));
    double *yx = malloc(m*sizeof(double));
    
    size_t nl = length(idx);

    size_t yp[nl+1];    
    size_t w0[nfeature];
    double x0[nfeature];
    memset(w0, 0, sizeof(size_t)*nfeature);
    memset(x0, 0, sizeof(double)*nfeature);

    size_t i;
    for (i = 0; i < nl; ++i) { // foreach cell
        yp[i] = n;
        size_t ii = INTEGER(idx)[i]-1;
        size_t j;
        for (j = wp[ii]; j < wp[ii+1]; ++j) {
            size_t k, p;
            size_t idx = wi[j];
            for (p = xp[idx]; p < xp[idx+1]; ++p) {
                k = xi[p];
                if (w0[k] < ii+1) {
                    if (n ==m) {
                        m = m*2;
                        Rprintf("%u\n", m);
                        yi = realloc(yi, sizeof(size_t)*m);
                        yx = realloc(yx, sizeof(double)*m);
                        if (yi == NULL || yx == NULL)
                            error("Failed to allocate data.");
                    }
                    
                    w0[k] = ii+1;
                    yi[n] = k;
                    x0[k] =wx[j]*xx[p];
                    n++;
                } else {
                    x0[k] += wx[j]*xx[p];
                }                                      
            }
            // n = scatter(x, wi[j], wx[j], w0, x0, ii+1, yi, n);
        }
        for (j = yp[i]; j < n; ++j) yx[j] = x0[yi[j]];
    }
    yp[i] = n;
    
    CHM_SP ans = M_cholmod_allocate_sparse(nfeature, nl, n, FALSE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) return(mkString("Failed to create sparse matrix"));
    size_t *ap = (size_t *)ans->p;
    size_t *ai = (size_t *)ans->i;
    double *ax = (double *)ans->x;
    memcpy(ap, yp, sizeof(size_t)*(nl+1));
    memcpy(ai, yi, sizeof(size_t)*n);
    memcpy(ax, yx, sizeof(double)*n);
    free(yi);
    free(yx);

    M_cholmod_sort(ans, &c);
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "N", R_NilValue);
}
