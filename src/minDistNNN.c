#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

double R_euclidean(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist;
    int count, j;
    count = 0;
    dist = 0;

    if (i1 >= nr) {
        error("i1 should smaller than nrow.");
    }

    if (i2 >= nr) {
        error("i2 should smaller than nrow.");
    }

    for (j = 0; j < nc; j++) {
        if (!ISNAN(x[i1]) && !ISNAN(x[i2])) {
            dev = (x[i1] - x[i2]);
            if (!ISNAN(dev)) {
                dist += dev * dev;
                count++;
            }
        }
        i1+=nr;
        i2+=nr;
    }
    if (count == 0) return NA_REAL;
    if (count != nc) dist /= ((double)count/nc);
    return sqrt(dist);
}
static int comp(const void* a, const void* b) {
    return (int)(*(double*)a - *(double*)b);
}
SEXP min_Dist_N_Nearest_Neighbors(SEXP x, SEXP _n)
{
    int n = asInteger(_n);
    int nr = nrows(x);
    int nc = ncols(x);
    
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
    
    // SEXP md = PROTECT(allocVector(REALSXP, nr));
    double min = 0;
    double d[nr];
    int i, j;
    for (i = 0; i < nr; ++i) {
        for (j = 0; j < nr; ++j) {
            if (i == j) d[j] = 0;
            else {
                d[j] = R_euclidean(REAL(x), nr, nc, i, j);
            }
        }
        qsort(d, nr, sizeof(double), comp);
        if (min == 0) min = d[n-1];
        if (min > d[n-1]) min = d[n-1];
        //REAL(md)[i] = d[n-1];
    }

    UNPROTECT(1);
    return ScalarReal(min);
}
