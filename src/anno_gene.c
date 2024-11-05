#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "gtf.h"

SEXP tx2gene(SEXP _tx, SEXP _db)
{
    int l = Rf_length(_tx);
    struct gtf_spec *G = (struct gtf_spec *)R_ExternalPtrAddr(_db);
    SEXP gene = PROTECT(allocVector(STRSXP, l));
    for (int i = 0; i < l; ++i) {
        const char *tx = translateChar(STRING_ELT(_tx, i));

        struct gtf *tt = gtf_query_tx(G, tx);
        if (tt == NULL) {
            SET_STRING_ELT(gene, i, mkChar("."));
        } else {
            char *gn = GTF_genename(G, tt->gene_name);
            SET_STRING_ELT(gene, i, mkChar(gn));
        }
    }
    UNPROTECT(1);
    return gene;
}


