#include<Rdefines.h>

SEXP build_weight_i(SEXP _n, SEXP _s)
{
    int n = asInteger(_n);
    int s = asInteger(_s);
    SEXP a = PROTECT(allocVector(INTSXP, n*s));
    int k = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < s; ++j) {
            int l = i +j - s/2;
            if (l < 0) l = 0;
            if (l >= n) l = n-1;
            INTEGER(a)[k++] = l;
        }
    }
    UNPROTECT(1);
    return a;
}

    
