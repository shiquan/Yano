#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "dict.h"
#include "gtf.h"
#include "bed.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"


struct val0 {
    union {
        double f;
        int64_t b;
    } a;
    char *c;
};

struct val {
    int id;
    int type;
    int number;
    int convert2str;
    struct val0 *v;
};
// edit from vcf.c
int process_fmt_array(kstring_t *s, int n, int type, void *data)
{
    int j = 0;
    uint32_t e = 0;
    if (n == 0) {
        return kputc('.', s) >= 0 ? 0 : -1;
    }
    if (type == BCF_BT_CHAR)
    {
        char *p = (char*)data;
        for (j = 0; j < n && *p; ++j, ++p)
        {
            if ( *p==bcf_str_missing ) e |= kputc('.', s) < 0;
            else e |= kputc(*p, s) < 0;
        }
    }
    else
    {
        #define BRANCH(type_t, convert, is_missing, is_vector_end, kprint) { \
            uint8_t *p = (uint8_t *) data; \
            for (j=0; j<n; j++, p += sizeof(type_t))    \
            { \
                type_t v = convert(p); \
                if ( is_vector_end ) break; \
                if ( j ) kputc(',', s); \
                if ( is_missing ) kputc('.', s); \
                else e |= kprint < 0; \
            } \
        }
        switch (type) {
            case BCF_BT_INT8:  BRANCH(int8_t,  le_to_i8, v==bcf_int8_missing,  v==bcf_int8_vector_end,  kputw(v, s)); break;
            case BCF_BT_INT16: BRANCH(int16_t, le_to_i16, v==bcf_int16_missing, v==bcf_int16_vector_end, kputw(v, s)); break;
            case BCF_BT_INT32: BRANCH(int32_t, le_to_i32, v==bcf_int32_missing, v==bcf_int32_vector_end, kputw(v, s)); break;
            case BCF_BT_FLOAT: BRANCH(uint32_t, le_to_u32, v==bcf_float_missing, v==bcf_float_vector_end, kputd(le_to_float(p), s)); break;
            default: hts_log_error("Unexpected type %d", type); exit(1); break;
        }
        #undef BRANCH
    }
    return e == 0 ? 0 : -1;
}

SEXP anno_vcf(SEXP _chr, SEXP _st, SEXP _ed, SEXP _ref, SEXP _alt, SEXP _strand, SEXP _vcf, SEXP _tags, SEXP check_alt_only)
{
    int l = Rf_length(_chr);
    if (Rf_length(_st) != l) {
        Rprintf("Inconsistance length of chr and start position.");
        return R_NilValue;
    }

    Rboolean check_alt = asLogical(check_alt_only);
    
    const char *vcf_fname = translateChar(STRING_ELT(_vcf, 0));

    htsFile *fp = hts_open(vcf_fname, "r");
    if (fp == NULL) {
        Rprintf("Failed to read %s.\n", vcf_fname);
        return R_NilValue;
    }
    
    htsFormat type0 = *hts_get_format(fp);    
    if ( type0.format != vcf && type0.format != bcf ) {
        Rprintf("Unsupport input format, only accept VCF/BCF.\n");
        return R_NilValue;
    }
    
    if ( fp->format.compression != bgzf ) {
        Rprintf("This file is NOT compressed by bgzip. %s", vcf_fname);
        return R_NilValue;
    }

    bcf_hdr_t *h = bcf_hdr_read(fp);    
    if (h == NULL) {
        Rprintf("Failed to load header of %s.\n", vcf_fname);
        return R_NilValue;
    }

    int n_tag = Rf_length(_tags);
    char **tags = calloc(n_tag, sizeof(char*));
    struct val *vals = calloc(n_tag, sizeof(struct val));

    kstring_t tmpk = {0,0,0};
    bcf1_t *v = bcf_init();

    hts_idx_t *idx = NULL;
    tbx_t *tbx = NULL;
    
    if (type0.format == bcf) {
        idx = bcf_index_load(vcf_fname);
        if (!idx) {
            Rprintf("Failed to load index file of %s. Use `bcftools index` the BCF first.\n", vcf_fname);
            return R_NilValue;
        }
    } else {
        tbx = tbx_index_load(vcf_fname);
        if (!tbx) {
            fprintf(stderr,"Failed to load index file of %s.\n", vcf_fname);
            return R_NilValue;
        }
    }
    
    for (int i = 0; i < n_tag; ++i) {
        tags[i] = strdup(translateChar(STRING_ELT(_tags, i)));
        
        int id = bcf_hdr_id2int(h, BCF_DT_ID, tags[i]);
        if (! bcf_hdr_idinfo_exists(h, BCF_HL_INFO, id) ) {
            Rprintf("No tag %s found in the VCF header.\n", tags[i]);
            return R_NilValue;
        }
        vals[i].id = id;
        vals[i].type = bcf_hdr_id2type(h, BCF_HL_INFO, id);
        vals[i].number = bcf_hdr_id2length(h, BCF_HL_INFO, id);
        vals[i].convert2str = 0;
        vals[i].v = calloc(l, sizeof(struct val0));
    }

    kstring_t str = {0,0,0};
    for (int i = 0; i < l; ++i) {
        const char *chr = translateChar(STRING_ELT(_chr, i));
        int start = INTEGER(_st)[i];
        
        const char *ref = translateChar(STRING_ELT(_ref, i));
        const char *alt = translateChar(STRING_ELT(_alt, i));
        
        for (int j = 0; j < n_tag; ++j) {
            struct val *val = &vals[j];
            val->v[i].a.f = bcf_float_missing;   
        }
        
        int rlen = strlen(ref);
        int tid;
        
        hts_itr_t *itr;
        if (type0.format == bcf) {
            tid = bcf_hdr_name2id(h, chr);
            itr = bcf_itr_queryi(idx, tid, start-1, start+rlen-1);
        } else {
            tid = tbx_name2id(tbx, chr);
            itr = tbx_itr_queryi(tbx, tid, start-1, start+rlen-1);
        }
        if (!itr) continue;
            
        for ( ;; ) {
            int ret;
            if (type0.format == bcf) {
                ret = bcf_itr_next(fp, itr, v);
                if (ret < 0) break;
            } else {
                str.l = 0;
                ret = tbx_itr_next(fp, tbx, itr, &str);
                if (ret < 0) break;
                vcf_parse(&str, h, v);
            }
            if (rlen != v->rlen) continue;
            if (start != v->pos+1) continue;

            bcf_unpack(v, BCF_UN_INFO);
            
            if (strcmp(v->d.allele[0], ref) != 0) {
                Rprintf("Inconsistant reference, make sure you use the right genome reference. tid: %d. %s:%d,%s vs %s, %d\n", tid, chr, start, ref, v->d.allele[0], v->rid);
                continue;
            }
            int allele = 0;
            if (check_alt) allele = 1;
            for (; allele < v->n_allele; allele++) {
                if (strcmp(v->d.allele[allele], alt) == 0) break;
            }
            
            if (allele == v->n_allele) continue;
            
            for (int j = 0; j < n_tag; ++j) {
                tmpk.l = 0;
                struct val *val = &vals[j];
                bcf_info_t *inf = bcf_get_info_id(v, val->id);
                
                if (inf==NULL) continue;
                if (allele == 0 && val->number ==  BCF_VL_A) continue;
                if (inf->len <= 0) {
                    val->v[i].a.b = 1;
                    continue;
                }
                
                if (inf->len == 1) {
                    switch (inf->type) {
                    case BCF_BT_INT8:
                        if (inf->v1.i != bcf_int8_missing ) val->v[i].a.b = inf->v1.i;
                        break;

                    case BCF_BT_INT16:
                        if (inf->v1.i != bcf_int16_missing ) val->v[i].a.b = inf->v1.i;
                        break;
                        
                    case BCF_BT_INT32:
                        if (inf->v1.i != bcf_int32_missing ) val->v[i].a.b = inf->v1.i;
                        break;

                    case BCF_BT_FLOAT:
                        if (!bcf_float_is_missing(inf->v1.f) ) val->v[i].a.f = inf->v1.f;
                        break;
                        
                    case BCF_BT_CHAR:
                        if (inf->v1.i != bcf_str_missing ) val->v[i].a.b = inf->v1.i;
                        break;
                        
                    default:
                        Rprintf("todo: type %d\n", inf->type);
                        break;
                    }
                } else {

                    if (val->number == BCF_VL_R || val->number == BCF_VL_A) {
                        int dst = -1;
                        if (val->number == BCF_VL_R) dst = allele;
                        else if (val->number == BCF_VL_A) dst = allele == -1 ? -1 : allele-1;

                        if (dst == -1) {
                            tmpk.l = 0;
                            process_fmt_array(&tmpk, inf->len, inf->type, inf->vptr);
                            val->v[i].c = strdup(tmpk.s);
                            val->convert2str = 1;
                        } else {
                            int j;
                            uint8_t *data = inf->vptr;
                            
                            if (inf->type == BCF_BT_CHAR) {
                                char *p = (char*)data;
                                for (j = 0; j < allele && *p; ++j, ++p);
                                if (j == allele) {
                                    tmpk.l = 0;
                                    kputc(*p, &tmpk);
                                    kputs("", &tmpk);
                                    val->v[i].c = strdup(tmpk.s);
                                }
                            } else {
                                
#define BRANCH(type_t, convert, is_vector_end, replace) {       \
                                    uint8_t *p = (uint8_t*)data;        \
                                    for (j = 0; j < allele; j++, p += sizeof(type_t)) { \
                                        type_t v0 = convert(p);         \
                                        if (is_vector_end) break;       \
                                    }                                   \
                                    if (j == allele) {                  \
                                        type_t v0 = convert(p);         \
                                        replace = v0;                   \
                                    }                                   \
                                }
                                
                                switch (inf->type) {
                                case BCF_BT_INT8:
                                    BRANCH(int8_t, le_to_i8, v0==bcf_int8_vector_end, val->v[i].a.b);
                                    break;
                                    
                                case BCF_BT_INT16:
                                    BRANCH(int16_t, le_to_i16, v0==bcf_int16_vector_end, val->v[i].a.b);
                                    break;
                                    
                                case BCF_BT_INT32:
                                    BRANCH(int32_t, le_to_i32, v0==bcf_int32_vector_end, val->v[i].a.b);
                                    break;
                                    
                                case BCF_BT_FLOAT:
                                    BRANCH(float, le_to_float, v0==bcf_float_vector_end, val->v[i].a.f);
                                    break;
                                default:
                                    Rprintf("Unexpected type %d.", inf->type);
                                    break;
                                }
#undef BRANCH
                            }
                        }
                    } else {
                        tmpk.l = 0;
                        process_fmt_array(&tmpk, inf->len, inf->type, inf->vptr);
                        val->v[i].c = strdup(tmpk.s);
                        val->convert2str = 1;
                    }
                }
            }
            break;
        }
        
        hts_itr_destroy(itr);
    }

    bcf_destroy(v);
    bcf_hdr_destroy(h);

    if (idx) hts_idx_destroy(idx);
    if (tbx) tbx_destroy(tbx);
    
    hts_close(fp);

    SEXP sl = PROTECT(allocVector(VECSXP, n_tag));
    for (int i = 0; i < n_tag; ++i) {
        struct val *val = &vals[i];
        if (val->type == BCF_HT_FLAG) {
            SEXP v = PROTECT(allocVector(INTSXP, l));
            for (int j = 0; j < l; ++j) {
                if (val->v[j].a.f == bcf_float_missing) {
                    INTEGER(v)[j] = NA_INTEGER;
                } else {
                    INTEGER(v)[j] = val->v[j].a.b;
                }
            }
            SET_VECTOR_ELT(sl, i, v);          
        } else if (val->type == BCF_HT_INT) {
            if (val->convert2str) {
                tmpk.l = 0;
                SEXP v = PROTECT(allocVector(STRSXP, l));
                for (int j = 0; j < l; ++j) {
                    if (val->v[j].c != NULL) {
                        SET_STRING_ELT(v, j, mkChar((const char*)val->v[j].c));
                        free(val->v[j].c);
                    }
                }
                SET_VECTOR_ELT(sl, i, v);
            } else {
                SEXP v = PROTECT(allocVector(INTSXP, l));
                for (int j = 0; j < l; ++j) {
                    if (val->v[j].a.f == bcf_float_missing) {
                        INTEGER(v)[j] = NA_INTEGER;
                    } else {
                        INTEGER(v)[j] = val->v[j].a.b;
                    }
                }
                SET_VECTOR_ELT(sl, i, v);
            }
        } else if (val->type == BCF_HT_REAL) {
            if (val->convert2str) {
                SEXP v = PROTECT(allocVector(STRSXP, l));
                for (int j = 0; j < l; ++j) {
                    if (val->v[j].c != NULL) {
                        SET_STRING_ELT(v, j, mkChar((const char*)val->v[j].c));
                        free(val->v[j].c);
                    }
                }
                SET_VECTOR_ELT(sl, i, v);
            } else {
                SEXP v = PROTECT(allocVector(REALSXP, l));
                for (int j = 0; j < l; ++j) {
                    if (val->v[j].a.f == bcf_float_missing) {
                        REAL(v)[j] = NA_REAL;
                    } else {
                        REAL(v)[j] = val->v[j].a.f;
                    }
                }
                SET_VECTOR_ELT(sl, i, v);
            }
        } else if (val->type == BCF_HT_STR) {
            SEXP v = PROTECT(allocVector(STRSXP, l));
            for (int j = 0; j < l; ++j) {
                if (val->v[j].c == NULL) {
                    tmpk.l = 0;
                    if (val->v[j].a.f != bcf_float_missing) {
                        kputc((char)val->v[j].a.b, &tmpk);
                        kputs("", &tmpk);
                        SET_STRING_ELT(v, j, mkChar(tmpk.s));
                    } 
                } else {                    
                    SET_STRING_ELT(v, j, mkChar((const char*)val->v[j].c));
                    free(val->v[j].c);
                }
            }
            SET_VECTOR_ELT(sl, i, v);
        }
        free(val->v);
    }
    free(vals);
    for (int i = 0; i < n_tag; ++i) free(tags[i]);
    free(tags);
    
    if (tmpk.m) free(tmpk.s);
    if (str.m) free(str.s);

    UNPROTECT(n_tag+1);

    return sl;
}

SEXP anno_gene(SEXP _chr, SEXP _st, SEXP _ed, SEXP _ref, SEXP _alt, SEXP _strand, SEXP _db)
{
    int l = Rf_length(_chr);
    if (Rf_length(_st) != l) {
        Rprintf("Inconsistance length of chr and start position.");
        return R_NilValue;
    }
    
    struct gtf_spec *G = (struct gtf_spec *)R_ExternalPtrAddr(_db);
    SEXP sl = PROTECT(allocVector(VECSXP, 2));
    SEXP gene = PROTECT(allocVector(STRSXP, l));
    SEXP type = PROTECT(allocVector(STRSXP, l));

    kstring_t tmp = {0,0,0};
    
    for (int i = 0; i < l; ++i) {
        const char *chr = translateChar(STRING_ELT(_chr, i));
        int start = INTEGER(_st)[i];

        // consequence
        //const char *ref = translateChar(STRING_ELT(_ref, i));
        //const char *alt = translateChar(STRING_ELT(_alt, i));
        
        int strand = -1;
        if (!Rf_isNull(_strand)) {
            const char *s = translateChar(STRING_ELT(_strand, i));
            if (s[0] == '+') strand = 0;
            else if (s[0] == '-') strand = 1;
        }

        int k;
        struct anno0 *a = anno_bed_core(chr, start-1, start, strand, G, &k, 1, 1000, 0, 1000, 1000);
        if (k == 0) {
            SET_STRING_ELT(gene, i, mkChar("."));
            SET_STRING_ELT(type, i, mkChar("intergenic"));
        } else if (k == 1) {
            SET_STRING_ELT(gene, i, mkChar(GTF_genename(G, a[0].g->gene_name)));
            SET_STRING_ELT(type, i, mkChar(bed_typename(a[0].type)));
            free(a);
        } else if (k > 1) {
            tmp.l = 0;
            for (int j = 0; j < k; ++j) {
                if (j) kputc(',', &tmp);
                kputs(GTF_genename(G, a[j].g->gene_name), &tmp);
            }
            kputs("", &tmp);
            SET_STRING_ELT(gene, i, mkChar((const char*)tmp.s));
            
            if (a[0].type > 9) {
                SET_STRING_ELT(type, i, mkChar(bed_typename(a[0].type)));
            } else {
                SET_STRING_ELT(type, i, mkChar("multigenes"));
            }
            free(a);
        }
    }

    if (tmp.m) free(tmp.s);
    SET_VECTOR_ELT(sl, 0, gene);
    SET_VECTOR_ELT(sl, 1, type);
    UNPROTECT(3);
    return sl;
}


