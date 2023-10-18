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
SEXP matrix_distance(SEXP _m, SEXP method)
{
    static const char *valid[] = {"dgCMatrix", ""};
    int ctype = R_check_class_etc(_m, valid);
    if (ctype < 0) return(mkString("Invalid class."));
    
    CHM_SP M = AS_CHM_SP__(_m);
    M = M_cholmod_transpose(M, (int)M->xtype, &c);

    int ncell = M->nrow;
    int nfeature = M->ncol;
    const int * mp = (int*)M->p;
    const int * mi = (int*)M->i;
    const double * mx = (double*)M->x;

    int nz = 0;
    for (int i = 1; i < ncell; ++i) nz+=i;

    CHM_SP ans = M_cholmod_allocate_sparse(ncell, ncell, nz, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    int *ap = (int*)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;
    
    double *a = (double*)R_Calloc(nfeature, double);
    double *b = (double*)R_Calloc(nfeature, double);

    method_func func = select_method(method);
    
    int k = 0;
    int i;
    for (i = 0; i < ncell; ++i) {
        memset(a, 0, sizeof(double)*nfeature);
        ap[i] = k;
        
        if (mp[i] == mp[i+1]) {
            for (int j = i; j < ncell; ++j) {
                ax[k] = 0;
                ai[k] = j;
                k++;
            }
            continue;
        }

        for (int l = mp[i]; l < mp[i+1]; ++l) {
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
                b[mi[l]] = mx[l];
            }

            ax[k] = func(a, b, nfeature);
            ai[k] = j;
            k++;
        }
    }
    
    ap[i] = k;

    assert(k == nz);
    
    R_Free(a);
    R_Free(b);

    M_cholmod_free_sparse(&M, &c);
    //UNPROTECT(2);
    return M_chm_sparse_to_SEXP(ans, 1, 0, 0, "", R_NilValue);;
}
