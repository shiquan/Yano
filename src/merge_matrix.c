#include <R.h>
#include <Rinternals.h>

#include "utils.h"
#include "dict.h"
#include "Matrix/Matrix.h"

SEXP merge_matrix(SEXP x)
{
    int N = length(x);
    struct dict *rns = dict_init();
    struct dict *cns = dict_init();
    int nzmax = 0;
    static const char *valid[] = {"dgCMatrix", ""};
    for (int i = 0; i < N; ++i) {
        SEXP a = VECTOR_ELT(x, i);
        int ivalid = R_check_class_etc(a, valid);
	if (ivalid < 0) {
            return(mkString("Only support dgCMatrix now!"));
        }
        
        SEXP pi = GET_SLOT(a, install("i"));
                            
        nzmax += length(pi);
        SEXP dimnames = GET_SLOT(a, install("Dimnames"));

        SEXP RN = VECTOR_ELT(dimnames,0);
        SEXP CN = VECTOR_ELT(dimnames,1);
        
        int lr = length(RN);
        for (int j = 0; j < lr; ++j) {
            const char *r = CHAR(STRING_ELT(RN,j));
            dict_push(rns, r);
        }
        

        int lc = length(CN);        
        for (int j = 0; j < lc; ++j) {
            const char *r = CHAR(STRING_ELT(CN,j));
            dict_push(cns, r);
        }
    }

    int nrow = dict_size(rns);
    int ncol = dict_size(cns);

    SEXP class = PROTECT(R_do_MAKE_CLASS("dgTMatrix"));
    SEXP O = PROTECT(R_do_new_object(class));
    SEXP dim = PROTECT(GET_SLOT(O, install("Dim")));
    INTEGER(dim)[0] = nrow;
    INTEGER(dim)[1] = ncol;

    SEXP RN = PROTECT(allocVector(STRSXP, nrow));
    for (int i = 0; i < nrow; ++i) {
        SET_STRING_ELT(RN, i, mkChar(dict_name(rns, i)));
    }

    SEXP CN = PROTECT(allocVector(STRSXP, ncol));
    for (int i = 0; i < ncol; ++i) {
        SET_STRING_ELT(CN, i, mkChar(dict_name(cns, i)));
    }

    SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dimnames, 0, RN);
    SET_VECTOR_ELT(dimnames, 1, CN);

    SET_SLOT(O, install("Dimnames"), dimnames);
    
    SEXP Oj = PROTECT(allocVector(INTSXP, nzmax));
    SEXP Oi = PROTECT(allocVector(INTSXP, nzmax));
    SEXP Ox = PROTECT(allocVector(REALSXP, nzmax));

    double *pOx = R_Calloc(nzmax, double);
    int *pOi = R_Calloc(nzmax, int);
    int *pOj = R_Calloc(nzmax, int);
    int nz = 0;
    
    for (int i = 0; i < N; ++i) {
        SEXP _a = VECTOR_ELT(x, i);
        SEXP dimnames = GET_SLOT(_a, install("Dimnames"));
        SEXP RN = VECTOR_ELT(dimnames,0);
        SEXP CN = VECTOR_ELT(dimnames,1);
        
        CHM_SP a = AS_CHM_SP__(_a);

        const int *ap = a->p;
        const int *ai = a->i;
        const double *ax = a->x;
        
        for (int j = 0; j < a->ncol; ++j) {
            const char *nm2 = CHAR(STRING_ELT(CN,j));
            int c1 = dict_query(cns, nm2);
            for (int p0 = ap[j]; p0 < ap[j+1]; p0++) {
                int i0 = ai[p0];                
                const char *nm1 = CHAR(STRING_ELT(RN,i0));
                i0 = dict_query(rns, nm1);
                //continue;
                pOi[nz] = i0;
                pOj[nz] = c1;
                pOx[nz] = ax[p0];
                nz++;
            }
        }
    }

    memcpy(REAL(Ox), pOx, nzmax*sizeof(double));
    memcpy(INTEGER(Oi), pOi, nzmax*sizeof(int));
    memcpy(INTEGER(Oj), pOj, nzmax*sizeof(int));

    SET_SLOT(O, install("i"), Oi);
    SET_SLOT(O, install("j"), Oj);
    SET_SLOT(O, install("x"), Ox);
    
    dict_destroy(rns);
    dict_destroy(cns);

    R_Free(pOi);
    R_Free(pOj);
    R_Free(pOx);
    
    UNPROTECT(9);
    
    return O;
}
