#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include "Matrix.h"
#include "cholmod.h"
#include <assert.h>

static double Jaccard(double *a, double *b, int n)
{
    int x = 0, y = 0;

    for (int i = 0; i < n; ++i) {
        if (a[i] == 0 && b[i] == 0) continue;
        if (a[i] != 0 && b[i] != 0) x++;
        y++;        
    }
    if (y == 0) return 0;
    return (double)x/y;
}

typedef double (*method_func)(double *, double *, int);

static method_func select_method(SEXP method)
{
    const char *m = CHAR(STRING_ELT(method,0));
    if (strcmp(m, "Jaccard") == 0) {
        return &Jaccard;
    } else {
        error("Unknown method, %s", m);
    }
}

static cholmod_common c;

// distance between columns in sparseMatrix
SEXP matrix_distance(SEXP _m, SEXP method, SEXP _threads)
{
    int n_thread = asInteger(_threads);
    
    static const char *valid[] = {"dgCMatrix", ""};
    int ctype = R_check_class_etc(_m, valid);
    if (ctype < 0) return(mkString("Invalid class."));
    
    CHM_SP M = AS_CHM_SP__(_m);

    int ncell = M->ncol;
    int nfeature = M->nrow;
    const int *mp = (int*)M->p;
    const int *mi = (int*)M->i;
    const double *mx = (double*)M->x;

    R_xlen_t nz;
    nz = (R_xlen_t)ncell*(ncell-1)/2;
    if (nz > INT_MAX) return(mkString("Too large size."));
    
    CHM_SP ans = M_cholmod_allocate_sparse(ncell, ncell, nz, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) return(mkString("Failed to create sparse matrix"));
    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;
    
    method_func func = select_method(method);
    int i;
    
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < ncell; ++i) {
        double *a = (double*)R_Calloc(nfeature, double);
        double *b = (double*)R_Calloc(nfeature, double);        
        memset(a, 0, sizeof(double)*nfeature);

        size_t k = i*(ncell-1) + i - ((i+1)*i)/2;
        ap[i] = k;
        
        if (mp[i] == mp[i+1]) {
            for (int j = i+1; j < ncell; ++j) {
                ax[k] = 0;
                ai[k] = j;
                k++;
            }
            continue;
        }
        
        for (int l = mp[i]; l < mp[i+1]; ++l) {
            if (ISNAN(mx[l])) continue;
            a[mi[l]] = mx[l];
        }
        
        for (int j = i+1; j < ncell; ++j) {
            if (mp[j] == mp[j+1]) {
                ax[k] = 0;
                ai[k] = j;
                k++;
                continue;
            }
            
            memset(b, 0, sizeof(double)*nfeature);
            for (int l = mp[j]; l < mp[j+1]; ++l) {
                if (ISNAN(mx[l])) continue;
                b[mi[l]] = mx[l];
            }
            
            ax[k] = func(a, b, nfeature);
            ai[k] = j;
            k++;
        }
        R_Free(a);
        R_Free(b);
    }
    
    ap[ncell] = nz;
    return M_chm_sparse_to_SEXP(ans, 1, -1, 0, "N", R_NilValue);;
}
