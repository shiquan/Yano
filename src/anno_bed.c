#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "dict.h"
#include "gtf.h"
#include "bed.h"
#include "htslib/kstring.h"

SEXP anno_bed(SEXP _chr, SEXP _st, SEXP _ed, SEXP _strand, SEXP _db, SEXP _prom, SEXP _up, SEXP _down, SEXP _at_up, SEXP _at_down)
{
    Rboolean promoter = asLogical(_prom);
    int upstream = asInteger(_up);
    int downstream = asInteger(_down);
    int at_up = asInteger(_at_up);
    int at_down = asInteger(_at_down);
    
    int l = Rf_length(_chr);
    if (Rf_length(_st) != l) {
        Rprintf("Inconsistance length of chr and start position.");
        return R_NilValue;
    }
    if (Rf_length(_ed) != l) {
        Rprintf("Inconsistance length of chr and start position.");
        return R_NilValue;
    }
    struct gtf_spec *G = (struct gtf_spec *)R_ExternalPtrAddr(_db);

    SEXP sl = PROTECT(allocVector(VECSXP, 2));
    SEXP gene = PROTECT(allocVector(STRSXP, l));
    SEXP type = PROTECT(allocVector(STRSXP, l));

    kstring_t str = {0, 0, 0};
    
    for (int i = 0; i < l; ++i) {
        const char *chr = translateChar(STRING_ELT(_chr, i));
        int start = INTEGER(_st)[i];
        int end = INTEGER(_ed)[i];
        int strand = BED_STRAND_UNK;

        str.l = 0;
        
        if (!Rf_isNull(_strand)) {
            const char *s = translateChar(STRING_ELT(_strand, i));
            if (s[0] == '+') strand = BED_STRAND_FWD;
            else if (s[0] == '-') strand = BED_STRAND_REV;
        }

        int k;
        struct anno0 *a = anno_bed_core(chr, start, end, strand, G, &k, promoter, downstream, upstream, at_down, at_up);

        if (k == 0) {
            SET_STRING_ELT(type, i, mkChar("intergenic"));
            SET_STRING_ELT(gene, i, mkChar("."));
        } else if (k == 1) {
            int type0 = a[0].type;
            SET_STRING_ELT(type, i, mkChar(bed_typename(type0)));
            if (a[0].g) {
                SET_STRING_ELT(gene, i, mkChar(GTF_genename(G, a[0].g->gene_name)));
            } else {
                SET_STRING_ELT(gene, i, mkChar("."));
            }
        } else if (k > 1) {
            struct gtf *g0 = a[0].g;
            int type0 = a[0].type;
            int j;
            for (j = 1; j < k; ++j) {
                struct gtf *g1 = a[j].g;
                if (g0->gene_name == g1->gene_name) continue;
                break;
            }

            if (j == k) {
                SET_STRING_ELT(type, i, mkChar(bed_typename(type0)));
                if (a[0].g) {
                    SET_STRING_ELT(gene, i, mkChar(GTF_genename(G, a[0].g->gene_name)));
                } else {
                    SET_STRING_ELT(gene, i, mkChar("."));
                }
            } else {
                if (a[0].type > 9) {
                    type0 = a[0].type;
                } else {
                    type0 = BAT_MULTIGENES;
                }
                for (j = 0; j < k; ++j) {
                    struct gtf *g = a[j].g;
                    if (g) {
                        if (j) kputc(',', &str);
                        kputs(GTF_genename(G, g->gene_name), &str);
                    }
                }
                SET_STRING_ELT(type, i, mkChar(bed_typename(type0)));
                if (str.l) SET_STRING_ELT(gene, i, mkChar(str.s));
                else {
                    SET_STRING_ELT(gene, i, mkChar("."));
                }
            }
        }
    }

    if (str.m) free(str.s);

    SET_VECTOR_ELT(sl, 0, gene);
    SET_VECTOR_ELT(sl, 1, type);
    UNPROTECT(3);
    return sl;
}


