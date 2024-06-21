#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>

#include <assert.h>
#include "dict.h"
#include "gtf.h"
#include "htslib/kstring.h"

static int cmpfunc1 (const void *_a, const void *_b)
{
    const struct gtf *a = *(const struct gtf**)_a;
    const struct gtf *b = *(const struct gtf**)_b;
    if (a->seqname != b->seqname) return (a->seqname > b->seqname) - (a->seqname < b->seqname);
    if (a->start != b->start) return (a->start > b->start) - (a->start < b->start);
    return (a->end > b->end) - (a->end < b->end) ;
}
// type == 1, gene
// type == 2, transcript
SEXP gene_tracks0(const char *chr, int start, int end, struct gtf_spec *G, struct dict *genes)
{
    struct region_itr *itr = gtf_query(G, chr, start, end);
    if (itr == NULL) return R_NilValue;

    struct gtf **pool = NULL;
    int n = 0;
    int m = 0;
    
    int i;
    for (i = 0; i < itr->n; ++i) {
        struct gtf *g = (struct gtf *)itr->rets[i];
        if (start > g->end) continue;
        if (end < g->start) continue;

        if (genes) {
            char *gene = GTF_genename(G, g->gene_name);
            int idx = dict_query(genes, gene);
            if (idx == -1) continue;
        }
        
        for (int j = 0; j < g->n_gtf; ++j) {
            if (m == n) {
                if (n == 0) {
                    m = 4;
                    pool = malloc(m*sizeof(void*));
                } else {
                    m = m *2;
                    pool = realloc(pool, m*sizeof(void *));
                }
            }
            pool[n++] = g->gtf[j];
        }        
    }

    if (n == 0) return R_NilValue;
    
    qsort(pool, n, sizeof(struct gtf**), cmpfunc1);

    int idx[n];
    int j = 0;
    int last_ends[n];
    memset(idx, 0, n*sizeof(int));
    memset(last_ends, 0, sizeof(int)*n);

    int n_exon = 0;
    for (i = 0; i < n; ++i) {
        struct gtf *g = pool[i];
        n_exon += g->n_gtf;
        idx[i] = 0;
        for (int k = 0; k < j; ++k) {
            if (g->start > last_ends[k]) {
                idx[i] = k+1;
                last_ends[k] = g->end;
                break;
            }
        }

        if (idx[i] == 0) {
            last_ends[j] = g->end;
            idx[i] = ++j;
        }
    }

    SEXP st = PROTECT(allocVector(INTSXP, n + n_exon));
    SEXP ed = PROTECT(allocVector(INTSXP, n + n_exon));
    SEXP str = PROTECT(allocVector(STRSXP, n + n_exon));
    SEXP type = PROTECT(allocVector(INTSXP, n + n_exon));
    SEXP gn = PROTECT(allocVector(STRSXP, n + n_exon));
    SEXP tx = PROTECT(allocVector(STRSXP, n + n_exon));
    SEXP id = PROTECT(allocVector(INTSXP, n + n_exon));

    int k = 0;
    for (i = 0; i < n; ++i) {
        struct gtf *g = pool[i];
        INTEGER(st)[k] = g->start;
        INTEGER(ed)[k] = g->end;
        char *gene = GTF_genename(G, g->gene_name);
        char *trans = GTF_transid(G, g->transcript_id);
        
        SET_STRING_ELT(gn, k, mkChar(gene));
        SET_STRING_ELT(tx, k, mkChar(trans));
        SET_STRING_ELT(str, k, mkChar(g->strand == 0 ? "+" : "-"));
        INTEGER(id)[k] = idx[i];
        INTEGER(type)[k] = 1;
        k++;
        for (int j = 0; j < g->n_gtf; ++j) {
            struct gtf *ex = g->gtf[j];
            INTEGER(st)[k] = ex->start;
            INTEGER(ed)[k] = ex->end;

            SET_STRING_ELT(str, k, mkChar(ex->strand == 0 ? "+" : "-"));
            SET_STRING_ELT(gn, k, mkChar(gene));
            SET_STRING_ELT(tx, k, mkChar(trans));
            INTEGER(id)[k] = idx[i];
            INTEGER(type)[k] = 2;
            k++;
        }
    }
    assert(k == n + n_exon);

    free(pool);
    region_itr_destroy(itr);
    
    SEXP cls, nam, sl;
    PROTECT(sl = allocVector(VECSXP, 7));
    SET_VECTOR_ELT(sl, 0, st);
    SET_VECTOR_ELT(sl, 1, ed);
    SET_VECTOR_ELT(sl, 2, str);
    SET_VECTOR_ELT(sl, 3, type);
    SET_VECTOR_ELT(sl, 4, id);
    SET_VECTOR_ELT(sl, 5, gn);
    SET_VECTOR_ELT(sl, 6, tx);
    
    PROTECT(cls = allocVector(STRSXP, 1)); 
    SET_STRING_ELT(cls, 0, mkChar("data.frame"));
    classgets(sl, cls);

    PROTECT(nam = allocVector(STRSXP, 7));
    SET_STRING_ELT(nam, 0, mkChar("start"));
    SET_STRING_ELT(nam, 1, mkChar("end"));
    SET_STRING_ELT(nam, 2, mkChar("strand"));
    SET_STRING_ELT(nam, 3, mkChar("type"));
    SET_STRING_ELT(nam, 4, mkChar("idx"));
    SET_STRING_ELT(nam, 5, mkChar("gene"));
    SET_STRING_ELT(nam, 6, mkChar("transcript"));
    namesgets(sl, nam);

    SEXP rownam;
    PROTECT(rownam = allocVector(STRSXP, n+n_exon));

    for (k = 0; k < n + n_exon; ++k) {
        char nm[11];
        snprintf(nm, sizeof(nm), "%10d", k+1);
        SET_STRING_ELT(rownam, k, mkChar(nm));
    }

    setAttrib(sl, R_RowNamesSymbol, rownam);
    UNPROTECT(11);

    return sl;
}

SEXP gene_tracks(SEXP chr, SEXP start, SEXP end, SEXP db, SEXP genes)
{
    //if (Rf_isNull(db)) return R_NilValue;
    const char *chr0 = translateChar(STRING_ELT(chr, 0));
    int start0 = asInteger(start);
    int end0 =  asInteger(end);
    
    struct gtf_spec *G = (struct gtf_spec *)R_ExternalPtrAddr(db);
    struct dict *genes0 = NULL;
    if (!Rf_isNull(genes)) {
        int l = length(genes);
        genes0 = dict_init();
        
        for (int i = 0; i < l; ++i) {
            const char *s = translateChar(STRING_ELT(genes, i));
            dict_push(genes0, s);
        }   
    }

    SEXP sl = gene_tracks0(chr0, start0, end0, G, genes0);

    if (genes0) dict_destroy(genes0);
    
    return sl;
}

SEXP G_query_gene(SEXP db, SEXP gene)
{
    struct gtf_spec *G = (struct gtf_spec *)R_ExternalPtrAddr(db);
    const char *gene0 = translateChar(STRING_ELT(gene, 0));
    struct gtf *g = gtf_query_gene(G, gene0);
    if (g == NULL) return R_NilValue;
    if (g->ext)
        Rprintf("More than one record found for gene %s, just pick the first record..", gene0);

    SEXP sl, col1, col2, col3, col4;
    PROTECT(sl = allocVector(VECSXP, 4));
    PROTECT(col1 = allocVector(STRSXP, 1));
    PROTECT(col2 = allocVector(INTSXP, 1));
    PROTECT(col3 = allocVector(INTSXP, 1));
    PROTECT(col4 = allocVector(STRSXP, 1));
    
    SET_STRING_ELT(col1, 0, mkChar(GTF_seqname(G, g->seqname)));
    INTEGER(col2)[0] = g->start;
    INTEGER(col3)[0] = g->end;
    SET_STRING_ELT(col4, 0, mkChar(g->strand == 0 ? "+" : "-"));

    SET_VECTOR_ELT(sl, 0, col1);
    SET_VECTOR_ELT(sl, 1, col2);
    SET_VECTOR_ELT(sl, 2, col3);
    SET_VECTOR_ELT(sl, 3, col4);
    
    UNPROTECT(5);
    return sl;
}
