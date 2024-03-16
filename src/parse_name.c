#include "utils.h"
#include "number.h"
#include "htslib/kstring.h"

SEXP parse_var_names(SEXP name)
{
    int l = Rf_length(name);

    SEXP sl = PROTECT(allocVector(VECSXP, 5));
    SEXP chr = PROTECT(allocVector(STRSXP, l));
    SEXP start = PROTECT(allocVector(INTSXP, l));
    SEXP ref = PROTECT(allocVector(STRSXP, l));
    SEXP alt = PROTECT(allocVector(STRSXP, l));
    SEXP strand = PROTECT(allocVector(STRSXP, l));
    
    kstring_t str = {0,0,0};
    kstring_t tmp = {0,0,0};
    for (int i = 0; i < l; ++i) {        
        const char *n = translateChar(STRING_ELT(name, i));
        str.l = 0;
        kputs(n, &str);
        tmp.l = 0;
        
        int j;
        // chr
        for (j = 0; j < str.l; ) {
            if (str.s[j] == ':') {
                SET_STRING_ELT(chr, i, mkChar(tmp.s));
                j++;
                break;
            }
            kputc(str.s[j], &tmp);
            j++;
        }

        if (j == str.l) {
            Rprintf("Not a EAT name, %s", str.s);
            free(str.s);
            if (tmp.m) free(tmp.s);

            UNPROTECT(6);
            return R_NilValue;
        }
        
        tmp.l = 0;
        for (; j < str.l;) {
            if (isdigit(str.s[j])) kputc(str.s[j], &tmp);
            else break;
            j++;
        }

        if (tmp.l == 0) {
            Rprintf("Not a EAT name, %s", str.s);
            free(str.s);
            if (tmp.m) free(tmp.s);

            UNPROTECT(6);
            return R_NilValue;            
        }

        int st = str2int(tmp.s);

        INTEGER(start)[i] = st;

        tmp.l = 0;
        for (;j < str.l;) {
            if (str.s[j] == '=') {
                if (tmp.l == 0) {
                    Rprintf("Not a EAT name, %s", str.s);
                    free(str.s);
                    if (tmp.m) free(tmp.s);
                    
                    UNPROTECT(6);
                    return R_NilValue;            
                }
                SET_STRING_ELT(ref, i, mkChar(tmp.s));
            } else if (str.s[j] == '>') {
                if (tmp.l == 0) {
                    Rprintf("Not a EAT name, %s", str.s);
                    free(str.s);
                    if (tmp.m) free(tmp.s);
                    
                    UNPROTECT(6);
                    return R_NilValue;            
                }
                SET_STRING_ELT(ref, i, mkChar(tmp.s));
                tmp.l = 0;
            } else if (str.s[j] == '/') {
                SET_STRING_ELT(alt, i, mkChar(tmp.s));
                j++;
                if (str.s[j] == '-') {
                    SET_STRING_ELT(strand, i, mkChar("-"));
                } else if (str.s[j] == '+') {
                    SET_STRING_ELT(strand, i, mkChar("+"));
                } else if (str.s[j] == '.') {
                    SET_STRING_ELT(strand, i, mkChar("."));
                } else {
                    Rprintf("Not a EAT name, %s", str.s);
                    free(str.s);
                    if (tmp.m) free(tmp.s);
                    
                    UNPROTECT(6);
                    return R_NilValue;
                }
                break;
            } else {
                kputc(str.s[j], &tmp);
            }
            j++;
        }
    }


    free(str.s);
    if (tmp.m) free(tmp.s);

    SET_VECTOR_ELT(sl, 0, chr);
    SET_VECTOR_ELT(sl, 1, start);
    SET_VECTOR_ELT(sl, 2, ref);
    SET_VECTOR_ELT(sl, 3, alt);
    SET_VECTOR_ELT(sl, 4, strand);

    UNPROTECT(6);
    return sl;
}
