#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "dict.h"
#include "gtf.h"

SEXP gtf_genes(struct gtf_spec *G)
{
    int i;
    int n = 0;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        n += ctg->n_gtf;
    }
    
    SEXP chr_          = PROTECT(allocVector(STRSXP, n));
    SEXP start_        = PROTECT(allocVector(INTSXP, n));
    SEXP end_          = PROTECT(allocVector(INTSXP, n));
    SEXP biotype_      = PROTECT(allocVector(STRSXP, n));
    SEXP strand_       = PROTECT(allocVector(STRSXP, n));
    SEXP geneid_       = PROTECT(allocVector(STRSXP, n));
    SEXP genename_     = PROTECT(allocVector(STRSXP, n));

    int idx = 0;
    for (i = 0; i < dict_size(G->name); i++) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        int j;
        for (j = 0; j < ctg->n_gtf; ++j) {
            struct gtf *gtf = ctg->gtf[j];
            SET_STRING_ELT(chr_,      idx, mkChar((const char *)dict_name(G->name, i)));
            INTEGER(start_)[idx] = gtf->start;
            INTEGER(end_)[idx] = gtf->end;
            int q = -1;
            if (gtf->attr != NULL) {
                q = dict_query(gtf->attr, "gene_type");
                if (q == -1) q = dict_query(gtf->attr, "gene_biotype");
                if (q == -1) q = dict_query(gtf->attr, "biotype");
            }
            char *biotype = q == -1 ? "Unknown" : dict_query_value(gtf->attr, q);
            SET_STRING_ELT(biotype_,  idx, mkChar((const char*)biotype));
            SET_STRING_ELT(strand_,   idx, mkChar(gtf->strand == 0 ? "+" : "-"));
            SET_STRING_ELT(geneid_,   idx, mkChar((const char*)dict_name(G->gene_id, gtf->gene_id)));
            SET_STRING_ELT(genename_, idx, mkChar((const char*)dict_name(G->gene_name, gtf->gene_name)));
            idx++;
        }
    }

    SEXP df_ = PROTECT(allocVector(VECSXP, 7));

    SET_VECTOR_ELT(df_,  0, chr_);
    SET_VECTOR_ELT(df_,  1, start_);
    SET_VECTOR_ELT(df_,  2, end_);
    SET_VECTOR_ELT(df_,  3, biotype_);
    SET_VECTOR_ELT(df_,  4, strand_);
    SET_VECTOR_ELT(df_,  5, geneid_);
    SET_VECTOR_ELT(df_,  6, genename_);

    SET_CLASS(df_, mkString("data.frame"));

    SEXP names = PROTECT(allocVector(STRSXP, 7));
    SET_STRING_ELT(names,  0, mkChar("chr"));
    SET_STRING_ELT(names,  1, mkChar("start"));
    SET_STRING_ELT(names,  2, mkChar("end"));
    SET_STRING_ELT(names,  3, mkChar("biotype"));
    SET_STRING_ELT(names,  4, mkChar("strand"));
    SET_STRING_ELT(names,  5, mkChar("geneid"));
    SET_STRING_ELT(names,  6, mkChar("genename"));
    setAttrib(df_, R_NamesSymbol, names);

    SEXP rownames = PROTECT(allocVector(INTSXP, 2));
    SET_INTEGER_ELT(rownames, 0, NA_INTEGER);
    SET_INTEGER_ELT(rownames, 1, -n);
    setAttrib(df_, R_RowNamesSymbol, rownames);
    
    UNPROTECT(10);
    return df_;
}
SEXP gtf_transcript(struct gtf_spec *G)
{
    int i, j;
    int n = 0;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        for (j = 0; j < ctg->n_gtf; ++j) {
            struct gtf *gtf = ctg->gtf[j];
            n+= gtf->n_gtf;
        }
    }
    
    SEXP chr_          = PROTECT(allocVector(STRSXP, n));
    SEXP start_        = PROTECT(allocVector(INTSXP, n));
    SEXP end_          = PROTECT(allocVector(INTSXP, n));
    SEXP biotype_      = PROTECT(allocVector(STRSXP, n));
    SEXP strand_       = PROTECT(allocVector(STRSXP, n));
    SEXP txid_         = PROTECT(allocVector(STRSXP, n));
    SEXP geneid_       = PROTECT(allocVector(STRSXP, n));
    SEXP genename_     = PROTECT(allocVector(STRSXP, n));

    int idx = 0;
    for (i = 0; i < dict_size(G->name); i++) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        for (j = 0; j < ctg->n_gtf; ++j) { // gene
            struct gtf *gtf0 = ctg->gtf[j];
            int k;
            for (k = 0; k < gtf0->n_gtf; ++k) {
                struct gtf *gtf = gtf0->gtf[k];
                SET_STRING_ELT(chr_,      idx, mkChar((const char *)dict_name(G->name, i)));
                INTEGER(start_)[idx]      = gtf->start;
                INTEGER(end_)[idx]        = gtf->end;
                int q = -1;
                if (gtf->attr != NULL) {
                    q = dict_query(gtf->attr, "transcript_type");
                    if (q == -1) q = dict_query(gtf->attr, "transcript_biotype");
                    if (q == -1) q = dict_query(gtf->attr, "biotype");
                }

                char *biotype = q == -1 ? "Unknown" : dict_query_value(gtf->attr, q);

                SET_STRING_ELT(biotype_,  idx, mkChar((const char*)biotype));
                SET_STRING_ELT(strand_,   idx, mkChar(gtf->strand == 0 ? "+" : "-"));
                SET_STRING_ELT(txid_,     idx, mkChar((const char*)dict_name(G->transcript_id, gtf->transcript_id)));
                SET_STRING_ELT(geneid_,   idx, mkChar((const char*)dict_name(G->gene_id, gtf->gene_id)));
                SET_STRING_ELT(genename_, idx, mkChar((const char*)dict_name(G->gene_name, gtf->gene_name)));
                idx++;
            }
        }
    }

    assert(n == idx);
    SEXP df_ = PROTECT(allocVector(VECSXP, 8));

    SET_VECTOR_ELT(df_,  0, chr_);
    SET_VECTOR_ELT(df_,  1, start_);
    SET_VECTOR_ELT(df_,  2, end_);
    SET_VECTOR_ELT(df_,  3, biotype_);
    SET_VECTOR_ELT(df_,  4, strand_);
    SET_VECTOR_ELT(df_,  5, txid_);
    SET_VECTOR_ELT(df_,  6, geneid_);
    SET_VECTOR_ELT(df_,  7, genename_);

    SET_CLASS(df_, mkString("data.frame"));

    SEXP names = PROTECT(allocVector(STRSXP, 8));
    SET_STRING_ELT(names,  0, mkChar("chr"));
    SET_STRING_ELT(names,  1, mkChar("start"));
    SET_STRING_ELT(names,  2, mkChar("end"));
    SET_STRING_ELT(names,  3, mkChar("biotype"));
    SET_STRING_ELT(names,  4, mkChar("strand"));
    SET_STRING_ELT(names,  5, mkChar("txid"));
    SET_STRING_ELT(names,  6, mkChar("geneid"));
    SET_STRING_ELT(names,  7, mkChar("genename"));
    setAttrib(df_, R_NamesSymbol, names);

    SEXP rownames = PROTECT(allocVector(INTSXP, 2));
    SET_INTEGER_ELT(rownames, 0, NA_INTEGER);
    SET_INTEGER_ELT(rownames, 1, -n);
    setAttrib(df_, R_RowNamesSymbol, rownames);
    
    UNPROTECT(11);
    
    return df_;
}
SEXP gtf_exon(struct gtf_spec *G)
{
    int i, j, k;
    int n = 0;
    for (i = 0; i < dict_size(G->name); ++i) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        for (j = 0; j < ctg->n_gtf; ++j) {
            struct gtf *gtf0 = ctg->gtf[j];
            for (k = 0; k < gtf0->n_gtf; ++k) {
                struct gtf *gtf = gtf0->gtf[k];
                n += gtf->n_gtf;
            }
        }
    }
    
    SEXP chr_          = PROTECT(allocVector(STRSXP, n));
    SEXP start_        = PROTECT(allocVector(INTSXP, n));
    SEXP end_          = PROTECT(allocVector(INTSXP, n));
    SEXP strand_       = PROTECT(allocVector(STRSXP, n));
    SEXP txid_         = PROTECT(allocVector(STRSXP, n));
    SEXP geneid_       = PROTECT(allocVector(STRSXP, n));
    SEXP genename_     = PROTECT(allocVector(STRSXP, n));

    int idx = 0;
    for (i = 0; i < dict_size(G->name); i++) {
        struct gtf_ctg *ctg = dict_query_value(G->name, i);
        for (j = 0; j < ctg->n_gtf; ++j) { // gene
            struct gtf *gtf0 = ctg->gtf[j];
            for (k = 0; k < gtf0->n_gtf; ++k) { // tx
                struct gtf *gtf1 = gtf0->gtf[k];
                int l;
                for (l = 0; l < gtf1->n_gtf; ++l) { // exon
                    struct gtf *gtf = gtf1->gtf[l];
                    SET_STRING_ELT(chr_,      idx, mkChar((const char *)dict_name(G->name, i)));
                    INTEGER(start_)[idx] = gtf->start;
                    INTEGER(end_)[idx] = gtf->end;
                    SET_STRING_ELT(strand_,   idx, mkChar(gtf->strand == 0 ? "+" : "-"));
                    SET_STRING_ELT(txid_,     idx, mkChar((const char*)dict_name(G->transcript_id, gtf->transcript_id)));
                    SET_STRING_ELT(geneid_,   idx, mkChar((const char*)dict_name(G->gene_id, gtf->gene_id)));
                    SET_STRING_ELT(genename_, idx, mkChar((const char*)dict_name(G->gene_name, gtf->gene_name)));
                    idx++;
                }
            }
        }
    }

    assert(n == idx);

    SEXP df_ = PROTECT(allocVector(VECSXP, 7));

    SET_VECTOR_ELT(df_,  0, chr_);
    SET_VECTOR_ELT(df_,  1, start_);
    SET_VECTOR_ELT(df_,  2, end_);
    SET_VECTOR_ELT(df_,  3, strand_);
    SET_VECTOR_ELT(df_,  4, txid_);
    SET_VECTOR_ELT(df_,  5, geneid_);
    SET_VECTOR_ELT(df_,  6, genename_);

    SET_CLASS(df_, mkString("data.frame"));

    SEXP names = PROTECT(allocVector(STRSXP, 7));
    SET_STRING_ELT(names,  0, mkChar("chr"));
    SET_STRING_ELT(names,  1, mkChar("start"));
    SET_STRING_ELT(names,  2, mkChar("end"));
    SET_STRING_ELT(names,  3, mkChar("strand"));
    SET_STRING_ELT(names,  4, mkChar("txid"));
    SET_STRING_ELT(names,  5, mkChar("geneid"));
    SET_STRING_ELT(names,  6, mkChar("genename"));
    setAttrib(df_, R_NamesSymbol, names);

    SEXP rownames = PROTECT(allocVector(INTSXP, 2));
    SET_INTEGER_ELT(rownames, 0, NA_INTEGER);
    SET_INTEGER_ELT(rownames, 1, -n);
    setAttrib(df_, R_RowNamesSymbol, rownames);
    
    UNPROTECT(10);
    return df_;
}

SEXP C_gtf2db(SEXP filename)
{
    if (!Rf_isString(filename)) error_return("Input is not String");
    const char *file = translateChar(STRING_ELT(filename, 0));

    fprintf(stderr, "%s\n", file);
    struct gtf_spec *G = gtf_read(file, 1);
    //free(_filename);

    SEXP genes = PROTECT(gtf_genes(G));
    SEXP trans = PROTECT(gtf_transcript(G));
    SEXP exons = PROTECT(gtf_exon(G));

    SEXP ls_ = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ls_, 0, genes);
    SET_VECTOR_ELT(ls_, 1, trans);
    SET_VECTOR_ELT(ls_, 2, exons);

    //SET_CLASS(df_, mkString("list"));

    gtf_destroy(G);
    Rprintf("Done!\n");
    UNPROTECT(4);
    return ls_;
}
