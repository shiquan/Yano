#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>

SEXP imputation1(SEXP _x, SEXP idx, SEXP _W, SEXP _filter)
{
    cholmod_common c;
    M_R_cholmod_start(&c);        
    double filter = asReal(_filter);
    CHM_SP x = AS_CHM_SP__(_x);
    CHM_SP W = AS_CHM_SP__(_W);
    if (W->stype) return mkString("W cannot be symmetric");

    const int * xp = (int*)x->p;
    const int * xi = (int*)x->i;
    const double * xx = (double*)x->x;

    const int * wp = (int*)W->p;
    const int * wi = (int*)W->i;
    const double * wx = (double*)W->x;

    int n_cells = W->nrow;
    int nfeature = x->nrow;
    
    size_t n = 0, m = 1000;
    int *yi = malloc(m*sizeof(int));
    double *yx = malloc(m*sizeof(double));
    
    int nl = length(idx);

    int yp[nl+1];    
    int w0[nfeature];
    double x0[nfeature];
    int i0[nfeature];
    memset(w0, 0, sizeof(int)*nfeature);
    memset(x0, 0, sizeof(double)*nfeature);
    memset(i0, 0, sizeof(int)*nfeature);
    int i;
    int n1;
    for (i = 0; i < nl; ++i) { // foreach cell
        yp[i] = n;
        
        int ii = INTEGER(idx)[i]-1;
        int j;
        int n1 = 0;
        
        for (j = wp[ii]; j < wp[ii+1]; ++j) {
            int k, p;
            int idx = wi[j];
            for (p = xp[idx]; p < xp[idx+1]; ++p) {
                k = xi[p];
                if (w0[k] < ii+1) {                    
                    w0[k] = ii+1;
                    i0[n1] = k;
                    x0[k] = wx[j]*xx[p];
                    n1++;
                } else {
                    x0[k] += wx[j]*xx[p];
                }                                      
            }
        }
        for (j = 0; j < n1; ++j) {
            int k = i0[j];
            if (x0[k] < filter) continue;
            if (n ==m) {
                m = m*2;
                yi = realloc(yi, sizeof(int)*m);
                yx = realloc(yx, sizeof(double)*m);
                if (yi == NULL || yx == NULL)
                    error("Failed to allocate data.");
            }

            yi[n] = k;
            yx[n] = x0[k];
            n++;
        }
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

CHM_SP imputation0(CHM_SP x, CHM_SP W, double filter)
{
    cholmod_common c;
    M_R_cholmod_start(&c);
    if (x->ncol != W->nrow)
        error("Inconsistance dims for multiplication.");
    
    const int * xp = (int*)x->p;
    const int * xi = (int*)x->i;
    const double * xx = (double*)x->x;

    const int * wp = (int*)W->p;
    const int * wi = (int*)W->i;
    const double * wx = (double*)W->x;

    int n_cell = W->nrow;
    int nfeature = x->nrow;
    
    size_t n = 0, m = 1000;
    int *yi = malloc(m*sizeof(int));
    double *yx = malloc(m*sizeof(double));
    
    int yp[n_cell+1];    
    int w0[nfeature];
    double x0[nfeature];
    int i0[nfeature];
    memset(w0, 0, sizeof(int)*nfeature);
    memset(x0, 0, sizeof(double)*nfeature);
    memset(i0, 0, sizeof(int)*nfeature);
    int i;
    int n1;
    for (i = 0; i < n_cell; ++i) { // foreach cell
        yp[i] = n;
        int j;
        int n1 = 0;
        
        for (j = wp[i]; j < wp[i+1]; ++j) {
            int k, p;
            int idx = wi[j];
            for (p = xp[idx]; p < xp[idx+1]; ++p) {
                k = xi[p];
                if (w0[k] < i+1) {                    
                    w0[k] = i+1;
                    i0[n1] = k;
                    x0[k] = wx[j]*xx[p];
                    n1++;
                } else {
                    x0[k] += wx[j]*xx[p];
                }                                      
            }
        }
        for (j = 0; j < n1; ++j) {
            int k = i0[j];
            if (x0[k] < filter) continue;
            if (n ==m) {
                m = m*2;
                yi = realloc(yi, sizeof(int)*m);
                yx = realloc(yx, sizeof(double)*m);
                if (yi == NULL || yx == NULL)
                    error("Failed to allocate data.");
            }

            yi[n] = k;
            yx[n] = x0[k];
            n++;
        }
    }
    yp[i] = n;
    
    CHM_SP ans = M_cholmod_allocate_sparse(nfeature, n_cell, n, FALSE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) error("Failed to create sparse matrix");
    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;
    memcpy(ap, yp, sizeof(int)*(n_cell+1));
    memcpy(ai, yi, sizeof(int)*n);
    memcpy(ax, yx, sizeof(double)*n);
    free(yi);
    free(yx);

    M_cholmod_sort(ans, &c);
    return ans;
}
