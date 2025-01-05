#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "gtf.h"

SEXP gene_locs(SEXP _gene, SEXP _db)
{
    int l = Rf_length(_gene);
    struct gtf_spec *G = (struct gtf_spec *)R_ExternalPtrAddr(_db);

    SEXP sl = PROTECT(allocVector(VECSXP, 4));
    SEXP chr = PROTECT(allocVector(STRSXP, l));
    SEXP start = PROTECT(allocVector(INTSXP, l));
    SEXP end = PROTECT(allocVector(INTSXP, l));
    SEXP strand = PROTECT(allocVector(STRSXP, l));
    
    for (int i = 0; i < l; ++i) {
        const char *gene = translateChar(STRING_ELT(_gene, i));

        struct gtf *tt = gtf_query_gene(G, gene);        
        if (tt == NULL) {
            SET_STRING_ELT(chr, i, NA_STRING);
            INTEGER(start)[i] = NA_INTEGER;
            INTEGER(end)[i] = NA_INTEGER;
            SET_STRING_ELT(strand, i, NA_STRING);
        } else {
            char *sn = GTF_seqname(G, tt->seqname);

            SET_STRING_ELT(chr, i, mkChar(sn));
            INTEGER(start)[i] = tt->start;
            INTEGER(end)[i] = tt->end;

            SET_STRING_ELT(strand, i, tt->strand == GTF_STRAND_FWD ? mkChar("+") : mkChar("-"));
        }
    }

    SET_VECTOR_ELT(sl, 0, chr);
    SET_VECTOR_ELT(sl, 1, start);
    SET_VECTOR_ELT(sl, 2, end);
    SET_VECTOR_ELT(sl, 3, strand);

    UNPROTECT(5);
    return sl;
}


