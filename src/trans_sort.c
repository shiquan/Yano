#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>

#include <assert.h>

SEXP trans_sort(SEXP starts, SEXP ends)
{
    int l = length(starts);
    assert(length(ends) == l);

    int last_ends[l];
    SEXP idx = PROTECT(allocVector(INTSXP, l));
    memset(last_ends, 0, sizeof(int)*l);

    int j = 0;
    for (int i = 0; i < l; ++i) {
        INTEGER(idx)[i] = 0;
        for (int k = 0; k < j; ++k) {
            if (INTEGER(starts)[i] > last_ends[k]) {
                INTEGER(idx)[i] = k+1;
                last_ends[k] = INTEGER(ends)[i];
                break;
            }
        }
        if (INTEGER(idx)[i] == 0) {
            last_ends[j] = INTEGER(ends)[i];
            INTEGER(idx)[i] = ++j;
        }
    }
    UNPROTECT(1);
    return idx;
}
