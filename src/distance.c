#ifdef _OPENMP
#include<omp.h>
#endif

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
    if (y == 0) return 1;
    return 1- (double)x/y;
}

static double euclidean(double *a, double *b, int n) {
    int i;
    double dev;
    double dist = 0;
    for (i = 0; i < n; ++i) {
        dev = a[i] - b[i];
        if (!ISNAN(dev)) {
            dist += dev *dev;
        }
    }
    return sqrt(dist);
}
typedef double (*method_func)(double *, double *, int);

static method_func select_method(SEXP method)
{
    const char *m = CHAR(STRING_ELT(method,0));
    if (strcmp(m, "Jaccard") == 0) {
        return &Jaccard;
    } else if (strcmp(m, "euclidean") == 0) {
        return &euclidean;
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

double *matrix_i(double *x, int i, int nc, int nr, double *a)
{
    int j;
    for (j = 0; j < nc; ++j) {
        a[j] = x[i];
        i+=nr;
    }
    return a;
}

// distance between positions, with prune.distance
SEXP matrix_distance2(SEXP x, SEXP method, SEXP _dist)
{
    method_func func = select_method(method);
    double prune_dist = asReal(_dist);
    
    switch(TYPEOF(x)) {
    case LGLSXP:
    case INTSXP:
    case REALSXP:
        break;
    default:
        error("'x' must be numeric");
    }

    if (TYPEOF(x) != REALSXP) x = coerceVector(x, REALSXP);
    PROTECT(x);

    int nr = nrows(x);
    int nc = ncols(x);

    int n = 0;
    int m = 10000;
    int *pp = R_Calloc(nr+1, int);
    int *ii = R_Calloc(m, int);
    double *xx = R_Calloc(m, double);
    
    double a[nc];
    double b[nc];
    int i, j;
    for (i = 0; i < nr; ++i) {
        pp[i] = n;
        matrix_i(REAL(x), i, nc, nr, a);
        for (j = i + 1; j < nr; ++j) {
            matrix_i(REAL(x), j, nc, nr, b);
            double d = func(a, b, nc);
            if (d <= prune_dist) {
                if (n == m) {
                    m = m *2;
                    ii = R_Realloc(ii, m, int);
                    xx = R_Realloc(xx, m, double);
                }
                ii[n] = j;
                xx[n] = d;
                n++;
            }
        }
    }

    pp[nr] = n;
    
    CHM_SP ans = M_cholmod_allocate_sparse(nr, nr, n, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) return(mkString("Failed to create sparse matrix"));
    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;
    memcpy(ap, pp, sizeof(int)*(nr+1));
    memcpy(ai, ii, sizeof(int)*n);
    memcpy(ax, xx, sizeof(double)*n);

    R_Free(pp);
    R_Free(ii);
    R_Free(xx);

    UNPROTECT(1);
    return M_chm_sparse_to_SEXP(ans, 1, -1, 0, "N", R_NilValue);;
}

// distance between positions, with prune.distance
// exp is a sparse matrix
// x is a dense matrix of spatial coordinates
// dist is prune distance to define nearby cells/spots/bins
// method: 1, mean; 2, distance
SEXP smooth_exp_with_large_weight_matrix(SEXP _exp, SEXP x, SEXP _dist, SEXP method)
{
    double prune_dist = asReal(_dist);
    if (prune_dist < 0) error("prune distance is < 0.");
    
    int norm = asInteger(method);
    
    switch(TYPEOF(x)) {
    case LGLSXP:
    case INTSXP:
    case REALSXP:
        break;
    default:
        error("'x' must be numeric");
    }

    if(TYPEOF(x) != REALSXP) x = coerceVector(x, REALSXP);
    PROTECT(x);

    CHM_SP exp0 = AS_CHM_SP__(_exp);

    CHM_SP exp = M_cholmod_transpose(exp0, (int)exp0->xtype, &c);
    if (c.status < CHOLMOD_OK) return mkString("Out of memory");

    int ncell = exp->nrow;
    int nfeature = exp->ncol;
    
    const int *mp = (int*)exp->p;
    const int *mi = (int*)exp->i;
    const double *mx = (double*)exp->x;

    int nr = nrows(x); // cells
    int nc = ncols(x); // dims

    if (nr != ncell) error("Inconsistance cell number.");
    
    double a[nc];
    double b[nc];
    double tmp[nr];
    
    int n = 0;
    int m = 10000;
    int *pp = R_Calloc(ncell+1, int);
    int *ii = R_Calloc(m, int);
    double *xx = R_Calloc(m, double);
    pp[0] = 0;
    int i, j;
    for (i = 0; i < nr; ++i) {
        memset(tmp, 0, sizeof(double)*nr);
        matrix_i(REAL(x), i, nc, nr, a);
        int cnt = 0;
        for (j = 0; j < nr; ++j) {
            if (i == j) continue;
            matrix_i(REAL(x), j, nc, nr, b);
            double d = euclidean(a, b, nc);
            if (d <= prune_dist) {
                cnt++;
                tmp[j] = d;
            }
        }
        if (norm == 1) {
            int k;
            for (k = 0; k < nr && cnt > 0; ++k) {
                if (tmp[k] > 0) tmp[k] = 1.0/cnt;
            }
        } else if (norm == 2) {
            double sum = 0;
            int k;
            for (k = 0; k < nr && cnt > 0; ++k) {
                if (tmp[k] > 0) {
                    tmp[k] = 1/tmp[k];
                    sum += tmp[k];
                }
            }
            for (k = 0; k < nr && sum > 0; ++k) {
                if (tmp[k] > 0) tmp[k] = tmp[k]/sum;
            }
        } else {
            error("Unknown normalised method.");
        }
        // Rprintf("nn, %d\n", cnt);
        // smooth data
        int i1, i2;
        for (i1 = 0; i1 < nfeature && cnt > 0; ++i1) {
            double sm = 0;
            for (i2 = mp[i1]; i2 < mp[i1+1]; ++i2) {
                int idx = mi[i2];
                if (tmp[idx] > 0) {
                    double x0 = mx[i2];
                    sm += x0 * tmp[idx];
                }
            }
            if (sm > 0) {
                if (n == m) {
                    m = m *2;
                    ii = R_Realloc(ii, m, int);
                    xx = R_Realloc(xx, m, double);
                }
                xx[n] = sm;
                ii[n] = i1;
                n++;
            }
            // pp[i1+1] = n;
        }
        pp[i+1] = n;
    }
    Rprintf("n, %d\n", n);
    M_R_cholmod_start(&c);
    CHM_SP ans = M_cholmod_allocate_sparse(nfeature, ncell, n, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    if (ans == NULL) return(mkString("Failed to create sparse matrix"));
    int *ap = (int *)ans->p;
    int *ai = (int *)ans->i;
    double *ax = (double *)ans->x;
    memcpy(ap, pp, sizeof(int)*(ncell+1));
    memcpy(ai, ii, sizeof(int)*n);
    memcpy(ax, xx, sizeof(double)*n);

    R_Free(pp);
    R_Free(ii);
    R_Free(xx);

    M_cholmod_free_sparse(&exp, &c);

    //CHM_SP X = M_cholmod_transpose(ans, (int)ans->xtype, &c);
    //if (c.status < CHOLMOD_OK) return mkString("Out of memory");
    //M_cholmod_free_sparse(&ans, &c);

    M_cholmod_finish(&c);
    UNPROTECT(1);
    return M_chm_sparse_to_SEXP(ans, 1, -1, 0, "N", R_NilValue);
}

