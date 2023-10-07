#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "utils.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/tbx.h"
#include "htslib/kstring.h"
#include "coverage.h"
#include "dict.h"

SEXP depth2matrix(SEXP _fname, SEXP _name, SEXP _start, SEXP _end, SEXP _strand,
                  SEXP _split, SEXP _tag, SEXP _umi_tag,
                  SEXP _mapq_thres, SEXP _cells, SEXP _n_cell, SEXP _group, SEXP _group_names)
{
    if (!Rf_isString(_fname)) error_return("Input is not String");
    
    const char *fname = translateChar(STRING_ELT(_fname, 0));
    const char *name = translateChar(STRING_ELT(_name,0)); // chromosome
    int start = asInteger(_start);
    int end = asInteger(_end);

    int ignore_strand = 0;
    int strand = asInteger(_strand);
    if (strand > 1) strand = -1;
    if (strand == -2) {
        ignore_strand = 1;
        strand = -1;
    }
    
    Rboolean split_by_bc = Rf_asLogical(_split);

    const char *tag = NULL;
    if (!Rf_StringBlank(_tag))
        tag = translateChar(STRING_ELT(_tag, 0));
    const char *umi_tag = NULL;
    if (!Rf_StringBlank(_umi_tag)) umi_tag = translateChar(STRING_ELT(_umi_tag, 0));

    int mapq_thres = asInteger(_mapq_thres);
    int n_cell = asInteger(_n_cell);
    
    // init bam
    htsFile *fp = hts_open(fname, "r");
    if (fp == NULL) {
        Rprintf("Failed to open file %s\n", fname);
        return R_NilValue;
    }

    hts_idx_t *idx = sam_index_load(fp, fname);
    if (idx == NULL) {
        Rprintf("Failed to load index file. Use `samtools index` the bam file first.\n");
        return R_NilValue;
    }
    
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == NULL) {
        Rprintf("Malformed bam file.\n");
        return R_NilValue;
    }

    // coverage
    int tid = bam_name2id(hdr, name);
    if (tid < 0) {
        Rprintf("Not found chromosome %s in bam file.\n", name);
        return R_NilValue;
    }

    struct dict *bc = dict_init();
    int fix_barcodes = 0;
    int *alias = NULL;
    int alias_tag = 0;
    if (n_cell > 0) {
        int n1 = Rf_length(_cells);
        int n2 = Rf_length(_group);
        if (n1 != n2 || n1 != n_cell) {
            error_return("Inconsistance length of cell vector.");
        }
        int i;
        for (i = 0; i < n_cell; i++) {
            const char *name = translateChar(STRING_ELT(_cells,i));
            int ret = dict_query(bc, name);
            if (ret != -1) {
                Rprintf("Duplicate cell name records? %s", name);
                return NULL;
            }
            dict_push(bc, name);
        }
        alias_tag = 1;
        fix_barcodes = 1;
        split_by_bc = TRUE;
        alias = INTEGER(_group);
        
        int n3 = Rf_length(_group_names);
        for (i = 0; i < n_cell; ++i) {
            alias[i] = alias[i]-1; // convert 1-based to 0-based
            if (alias[i] > n3) error_return("Inconsistance alias names.");
        }
    }

    struct depth *d = bam2depth(idx, tid, start, end, strand, fp, mapq_thres, ignore_strand,
                                bc, tag, umi_tag, split_by_bc, alias_tag,
                                alias, fix_barcodes);
    // return matrix
    hts_close(fp);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);

    int n = 0;
    struct depth *r = d;
    
    for (; r; ) {
        r = r->next;
        n++;
    }

    if (n == 0) return R_NilValue;
    if (strand == -1) n = n*2;
    
    SEXP pos = PROTECT(allocVector(INTSXP, n));
    SEXP stra = PROTECT(allocVector(INTSXP, n));
    SEXP depth = PROTECT(allocVector(INTSXP, n));
    SEXP label = PROTECT(allocVector(STRSXP, n));

    int i = 0;
    r = d;
    for (; r; ) {
        INTEGER(pos)[i] = r->pos;
        INTEGER(stra)[i] = 0;
        INTEGER(depth)[i] = r->dep1;
        if (r->id == -1) {
            SET_STRING_ELT(label, i, mkChar("."));
        } else {
            if (alias_tag)
                SET_STRING_ELT(label, i, mkChar(translateChar(STRING_ELT(_group_names,r->id))));
            else
                SET_STRING_ELT(label, i, mkChar((const char*)dict_name(bc, r->id)));
        }

        i++;
        
        if (strand == -1) {
            INTEGER(pos)[i] = r->pos;
            INTEGER(stra)[i] = 1;
            INTEGER(depth)[i] = r->dep2;
            if (r->id == -1) {
                SET_STRING_ELT(label, i, mkChar("."));
            } else {
                if (alias_tag)
                    SET_STRING_ELT(label, i, mkChar(translateChar(STRING_ELT(_group_names,r->id))));
                else
                    SET_STRING_ELT(label, i, mkChar((const char*)dict_name(bc, r->id)));
            }
            i++;
        }
        
        struct depth *t = r;
        r = r->next;
        free(t);
    }
    assert(i == n);

    dict_destroy(bc);
    
    SEXP df = PROTECT(Rf_allocVector(VECSXP,4));
    SET_VECTOR_ELT(df, 0, pos);
    SET_VECTOR_ELT(df, 1, stra);
    SET_VECTOR_ELT(df, 2, depth);
    SET_VECTOR_ELT(df, 3, label);

    UNPROTECT(5);
    return df;
}

SEXP fragment2matrix(SEXP _fname, SEXP _name, SEXP _start, SEXP _end,
                     // SEXP _split,
                     SEXP _cells, SEXP _n_cell, SEXP _group, SEXP _group_names)
{
    if (!Rf_isString(_fname)) {
        Rprintf("Input is not String");
        return R_NilValue;
    }
    
    const char *fname = translateChar(STRING_ELT(_fname, 0));
    const char *seqname = translateChar(STRING_ELT(_name,0)); // chromosome
    int start = asInteger(_start);
    int end = asInteger(_end);
    
    //Rboolean split_by_bc = Rf_asLogical(_split);

    int n_cell = asInteger(_n_cell);

    BGZF *fp = bgzf_open(fname, "r");
    if (fp == NULL) {
        Rprintf("%s : %s.", fname, strerror(errno));
        return R_NilValue;
    }
            
    tbx_t *tbx = tbx_index_load(fname);
    if (tbx == NULL) {
        Rprintf("Failed to load index file.");
        return R_NilValue;
    }
    
    struct dict *bc = dict_init();
    // int fix_barcodes = 0;
    int *alias = NULL;
    // int alias_tag = 0;
    if (n_cell > 0) {
        int n1 = Rf_length(_cells);
        int n2 = Rf_length(_group);
        if (n1 != n2 || n1 != n_cell) {
            Rprintf("Inconsistance length of cell vector.");
            return R_NilValue;
        }

        for (int i = 0; i < n_cell; i++) {
            const char *name = translateChar(STRING_ELT(_cells,i));
            int ret = dict_query(bc, name);
            if (ret != -1) {
                Rprintf("Duplicate cell name records? %s", name);
                return R_NilValue;
            }
            dict_push(bc, name);
        }
        // alias_tag = 1;
        //fix_barcodes = 1;
        //split_by_bc = TRUE;
        alias = INTEGER(_group);
        
        int n3 = Rf_length(_group_names);
        for (int i = 0; i < n_cell; ++i) {
            alias[i] = alias[i]-1; // convert 1-based to 0-based
            if (alias[i] > n3) {
                Rprintf("Inconsistance alias names.");
                return R_NilValue;
            }
        }
    }

    struct depth *d = fragment2depth(tbx, seqname, start, end, fp, bc, alias);

    bgzf_close(fp);
    tbx_destroy(tbx); // ??

    int n = 0;
    struct depth *r = d;
    
    for (; r; ) {
        r = r->next;
        n++;
    }

    
    SEXP pos = PROTECT(allocVector(INTSXP, n));
    SEXP stra = PROTECT(allocVector(INTSXP, n));
    SEXP depth = PROTECT(allocVector(INTSXP, n));
    SEXP label = PROTECT(allocVector(STRSXP, n));

    int i = 0;
    r = d;
    for (; r; ) {
        INTEGER(pos)[i] = r->pos;
        INTEGER(stra)[i] = 0;
        INTEGER(depth)[i] = r->dep1;
        if (r->id == -1) {
            SET_STRING_ELT(label, i, mkChar("."));
        } else {
            if (alias)
                SET_STRING_ELT(label, i, mkChar(translateChar(STRING_ELT(_group_names,r->id))));
            else
                SET_STRING_ELT(label, i, mkChar((const char*)dict_name(bc, r->id)));
        }

        i++;
        struct depth *t = r;
        r = r->next;
        free(t);
    }
    assert(i == n);

    dict_destroy(bc);
    
    SEXP df = PROTECT(Rf_allocVector(VECSXP,4));
    SET_VECTOR_ELT(df, 0, pos);
    SET_VECTOR_ELT(df, 1, stra);
    SET_VECTOR_ELT(df, 2, depth);
    SET_VECTOR_ELT(df, 3, label);

    UNPROTECT(5);
    return df;
}

