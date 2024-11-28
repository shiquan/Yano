#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "dict.h"
#include "gtf.h"
#include "bed.h"
#include "variant_type.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"

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
// this function edit from bcftools/vcf.c
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
        Rprintf("Inconsistance length of chr and start position.\n");
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
                    SET_STRING_ELT(v, j, NA_STRING);
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
                    SET_STRING_ELT(v, j, NA_STRING);
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
                SET_STRING_ELT(v, j, NA_STRING);
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
        int strand = -1;
        if (!Rf_isNull(_strand)) {
            const char *s = translateChar(STRING_ELT(_strand, i));
            if (s[0] == '+') strand = BED_STRAND_FWD;
            else if (s[0] == '-') strand = BED_STRAND_REV;
        }

        int k;
        struct anno0 *a = anno_bed_core(chr, start-1, start, strand, G, &k, 1000, 1000);
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
                SET_STRING_ELT(type, i, mkChar(bed_typename(BAT_MULTIGENES)));
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

// Sequence ontology variant types
enum mol_con {
    mc_reference = 0,              // reference allele
    mc_whole_gene,                 // delete of whole gene
    mc_exon_loss,                  // delete of an exon
    mc_splice_donor,               // start region of an intron (2bp at the 5' end of intron)
    mc_splice_acceptor,            // end region of an intron (2bp at the 3' end of intron)
    mc_stop_gained,                // nonsense
    mc_exon_splice_sites,
    mc_frameshift_truncation,     
    mc_frameshift_elongation,      // translational reading frame to be extended relative to the reference
    mc_stop_loss,
    mc_start_loss,                 //
    mc_inframe_indel,
    mc_missense,                   //
    mc_synonymous,                 //
    mc_nocall,                     // variant different from reference genome but same with RNA reference
    mc_stop_retained,              //
    mc_start_retained,             //
    mc_utr5_premature_start_codon_gain,
    mc_utr5_exon,                  //
    mc_utr3_exon,                  //
    mc_splice_region,               // change with in 3-8 bases of the intron
    mc_coding_intron,
    mc_utr3_intron,
    mc_utr5_intron,    
    mc_noncoding_exon,
    mc_noncoding_splice_region,
    mc_noncoding_intron,
    mc_antisense_utr3,              // variant located at antisense UTR3 region
    mc_antisense_utr5,
    mc_antisense_exon,
    mc_antisense_intron,
    mc_upstream_1K,
    mc_downstream_1K,
    mc_antisense_upstream_1K,
    mc_antisense_downstream_1K,
    mc_intergenic,
    mc_unknown
};

const char *MCT[] = {
    "ref",                        //
    "whole_gene_delete",           //
    "exon_loss",                  //
    "splice_donor",               //
    "splice_acceptor",            //
    "stop_gained",                //
    "splice_sites",               //      
    "frameshift_truncation",      //
    "frameshift_elongation",      //
    "stop_loss",                  //
    "start_loss",                 //
    "inframe_indel",              //
    "missense",                    //
    "synonymous",                  //
    "nocall",                     // not implement yet
    "stop_retained",               //
    "start_retained",              //
    "utr5_start_gain",   // not implement yet
    "utr5",                        //
    "utr3",                        //
    "splice_region",               //
    "intron",                      //
    "utr3_intron",                 //
    "utr5_intron",                 //
    "noncoding_exon",              //
    "noncoding_splice_region",     //
    "noncoding_intron",            //
    "antisense_utr3",              //
    "antisense_utr5",              //
    "antisense_exon",              //
    "antisense_intron",            //
    "upstream_1K",                 //
    "downstream_1K",               //
    "antisense_upstream_1K",       //
    "antisense_downstream_1K",     //
    "intergenic",                   //
    "unknown"
};

struct csq {
    enum mol_con csq;
    int gene_name;
    struct gtf *exon;
};

int cmp(const void* a, const void* b) {
    struct csq *a1 = (struct csq*)a;
    struct csq *b1 = (struct csq*)b;
    return a1->csq - b1->csq;
}
struct mempool {
    struct dict *tx;
};

static struct mempool mempool;

void memory_init()
{
    mempool.tx = dict_init();
    dict_set_value(mempool.tx);
}
void memory_release()
{
    int i;
    for (i = 0; i < dict_size(mempool.tx); ++i) {
        void *v = dict_query_value(mempool.tx, i);
        free(v);
    }
    dict_destroy(mempool.tx);
}

#include "variant_type.h"

char *construct_reference(struct gtf *tx, struct gtf_spec *G, const char *fasta, faidx_t *fai)
{
    char *nm = GTF_transid(G, tx->transcript_id);
    char *v = dict_query_value2(mempool.tx, nm);
    if (v) return v;

    if (dict_size(mempool.tx) >= 1000) {
        memory_release();
        memory_init();
    }
    // Rprintf("cstart, %d, cend, %d\n", tx->cstart, tx->cend);
    int i;
    kstring_t tmp = {0,0,0};
    for (i = 0; i < tx->n_gtf; ++i) {
        struct gtf *ex = tx->gtf[i];
        if (ex->type != feature_exon) continue;
        int st = ex->start;
        int ed = ex->end;
        if (st > tx->cend) break;
        if (st < tx->cstart) st = tx->cstart;
        if (ed > tx->cend) ed = tx->cend;
        if (st < ed) {
            // Rprintf("start, %d, end, %d\n", st, ed);
            int len;
            char *seq = faidx_fetch_seq(fai, GTF_seqname(G, tx->seqname), st-1, ed-1, &len);
            if (seq == NULL) {
                if (tmp.m) free(tmp.s);
                return NULL;
            }
            kputs(seq, &tmp);
        }
    }

    if (tmp.l %3) {
        Rprintf("CDS region of transcript %s cannot translate to codons properly.\n",nm);
    }

    if (tx->strand == GTF_STRAND_REV) {
        compl_seq(tmp.s, tmp.l);        
    }
    int id = dict_push(mempool.tx, nm);
    dict_assign_value(mempool.tx, id, tmp.s);
    return tmp.s;
}
int find_cds_location(struct gtf *tx, int pos, int *ss)
{
    *ss = 0;
    if (pos < tx->cstart) return -1;
    if (pos > tx->cend) return -1;

    int i;
    int cds = 0;
    if (tx->strand == GTF_STRAND_REV) {
        for (i = tx->n_gtf-1; i >= 0; i--) {
            struct gtf *ex = tx->gtf[i];
            if (ex->type != feature_exon) continue;
            if (ex->start > tx->cend) continue;
            if (ex->end > tx->cend) {
                if (pos >= ex->start) {
                    if (pos-ex->start<3) *ss = 1;
                    return tx->cend - pos +1;   
                }
                cds = tx->cend - ex->start+1;
            } else {
                if (pos < ex->start) cds += ex->end - ex->start +1;
                else {
                    if (ex->end -pos < 3 || pos - ex->start<3) *ss = 1;
                    return ex->end - pos +1 + cds;
                }
            }
        }
    } else {
        for (i = 0; i < tx->n_gtf; ++i) {
            struct gtf *ex = tx->gtf[i];
            if (ex->type != feature_exon) continue;
            if (ex->end < tx->cstart) continue;
            if (ex->start < tx->cstart) {
                if (pos < ex->end) {
                    if (ex->end -pos <3) *ss = 1;
                    return pos - tx->cstart+1;   
                }
                cds += ex->end - tx->cstart+1;
            } else {
                if (pos > ex->end) cds += ex->end - ex->start+1;
                else {
                    if (ex->end -pos < 3 || pos -ex->start<3) *ss = 1;
                    return pos-ex->start+1 +cds;   
                }
            }
        }
    }
    return -1;
}
char *construct_alternative_sequence(char *ref_str, struct gtf *tx, int pos, const char *ref, const char *alt)
{
    int ss;
    int cds = find_cds_location(tx, pos, &ss);
    assert(cds>0);
    kstring_t tmp = {0,0,0};
    int lref = strlen(ref);
    int lalt = strlen(alt);

    if (tx->strand == GTF_STRAND_FWD) {        
        int i;
        for (i = 0; i < lref; ++i) {
            if (ref_str[cds-1+i] != ref[i]) {
                // Rprintf("Inconsist of reference allele and reference. Make sure you use the right fasta. %d, %s > %s\n", pos, ref, alt);
                // splice sits
                return NULL;
            }
        }
    } else {
        int i;
        for (i = 0; i < lref; ++i) {
            if (ref_str[cds-1+i] != revseqarr[seq2code4(ref[lref-1-i])]) {
                // Rprintf("Inconsist of reference allele and reference. Make sure you use the right fasta. %d, %s > %s\n", pos, ref, alt);
                // splice sits
                return NULL;
            }
        }

    }
    kputsn(ref_str, cds, &tmp);
    if (lalt > 0) {
        if (tx->strand == GTF_STRAND_FWD) kputs(alt, &tmp);
        else {
            int i;
            for (i = 0; i < lalt; ++i) {
                kputc(revseqarr[seq2code4(alt[lalt-i-1])], &tmp);
            }
            kputs("", &tmp);
        }
    }
    kputs(ref_str+cds+lref, &tmp);
        
    return tmp.s;
}
enum mol_con predict_molecular_consequence(char *ref, char *alt, int mt)
{
    int lref = strlen(ref);
    int lalt = strlen(alt);
    int i;
    for (i = 0; i < lref && i < lalt; i = i+3) {
        int ref_aa = codon2aminoid(ref+i, mt);
        int alt_aa = codon2aminoid(alt+i, mt);
        if (ref_aa != alt_aa) {
            if (ref_aa == C4_Stop) return mc_frameshift_elongation;
            if (alt_aa == C4_Stop) return mc_frameshift_truncation;
        }
    }
    
    if (lref == lalt) return mc_inframe_indel;
    return mc_frameshift_elongation; // unknown stop position, consider as elongation
}
enum mol_con fast_predict(int pos, struct gtf *tx, int *e)
{
    int i;
    int last_ed = -1;
    int last_i = -1;
    for (i = 0; i < tx->n_gtf; ++i) {
        struct gtf *ex = tx->gtf[i];
        if (ex->type != feature_exon) continue;
        if (pos >= ex->start && pos <= ex->end) return mc_noncoding_exon;
        if (pos < ex->start) {
            assert(last_i > 0);
            int loc1 = pos - last_ed;
            int loc2 = ex->start - pos;
            if (loc1 < loc2) {
                *e = last_i;
                if (loc1 < 9) return mc_noncoding_splice_region;
                return mc_noncoding_intron;
            } else {
                *e = i;
                if (loc2 < 9) return mc_noncoding_splice_region;
                return mc_noncoding_intron;
            }
        }
        if (pos > ex->end) {
            last_ed = ex->end;
            last_i = i;
        }
    }
    return mc_unknown;
}

enum mol_con predict_func0(const char *chr, int pos, int strand, const char *ref, const char *alt, struct gtf_spec *G, struct gtf *tx, const char *fasta, faidx_t *fai, int debug)
{
    enum mol_con ret = mc_reference;
    int lref = strlen(ref);
    int lalt = strlen(alt);
    if (lref == lalt && strcmp(ref, alt) == 0) return ret;
    int e = 0;
    if (tx->cstart == tx->cend) {
        ret = fast_predict(pos,tx, &e);
        if (strand != -1 && strand != tx->strand) {
            if (ret == mc_noncoding_exon) return mc_antisense_exon;
            return mc_antisense_intron;
        }
    } else {
        ret = fast_predict(pos, tx, &e);
        if (strand != -1 && strand != tx->strand) {
            if (ret == mc_noncoding_exon) {
                if (pos < tx->start)
                    return tx->strand == GTF_STRAND_FWD ? mc_antisense_utr5 : mc_antisense_utr3;
                if (pos > tx->end)
                    return tx->strand == GTF_STRAND_FWD ? mc_antisense_utr3 : mc_antisense_utr5;
                return mc_antisense_exon;
            }
            return mc_antisense_intron;
        }
        struct gtf *ex = tx->gtf[e];
        if (lref - lalt >= 10 && lref-lalt >= ex->end - ex->start) return mc_exon_loss;
        
        if (ret == mc_noncoding_intron) {
            if (pos < tx->start)
                return tx->strand == GTF_STRAND_FWD ? mc_utr5_intron : mc_utr3_intron;
            if (pos > tx->end)
                return tx->strand == GTF_STRAND_FWD ? mc_utr3_intron : mc_utr5_intron;
            return mc_coding_intron;            
        }
        if (ret == mc_noncoding_splice_region) {
            if (pos < ex->start && ex->start - pos < 3)
                return tx->strand == GTF_STRAND_FWD ? mc_splice_donor : mc_splice_acceptor;
            if (pos > ex->end && pos - ex->end < 3)
                return tx->strand == GTF_STRAND_FWD ? mc_splice_donor : mc_splice_acceptor;
            return mc_splice_region;
        }

        if (pos < tx->cstart) {
            return tx->strand == GTF_STRAND_FWD ? mc_utr5_exon : mc_utr3_exon;
        }
        if (pos > tx->cend) {
            return tx->strand == GTF_STRAND_FWD ? mc_utr3_exon : mc_utr5_exon;
        }
        
        char *ref_str = construct_reference(tx, G, fasta, fai);
        if (ref_str == NULL ) {
            return mc_exon_splice_sites;
            // continue; // return mc_unknown;   
        }
        
        if (lref == 1 && lalt == 1) {
            int ss;
            int cds = find_cds_location(tx, pos, &ss);
            char codon[3];
            codon[0] = *(ref_str + (cds-1)/3*3);
            codon[1] = *(ref_str + (cds-1)/3*3+1);
            codon[2] = *(ref_str + (cds-1)/3*3 +2);            
            int ref_aa = codon2aminoid(codon, tx->is_mito);

            if (debug) {
                Rprintf("transcript: %s\nCDS : %d\ncodon : %s\n", GTF_transid(G, tx->transcript_id), cds, codon);
            }
            if (tx->strand == GTF_STRAND_FWD) {
                codon[(cds-1) %3] = *alt;
            } else {
                codon[(cds-1) %3] = revseqarr[seq2code4(*alt)];
            }
            int alt_aa = codon2aminoid(codon, tx->is_mito);
            if (debug) {
                Rprintf("mut : %s\n%s > %s\n", codon, codon2name(ref_aa), codon2name(alt_aa));
            }
            if (ref_aa == alt_aa) {
                if (cds < 4) return mc_start_retained;
                if (ref_aa == C4_Stop) return mc_stop_retained;
                return mc_synonymous;
            } else {
                if (cds < 4) return mc_start_loss;
                if (alt_aa == C4_Stop) return mc_stop_gained;
                if (ss) return mc_exon_splice_sites;
                if (ref_aa == C4_Stop) return mc_stop_loss;
                return mc_missense;
            }

        } else {
            char *alt_str = construct_alternative_sequence(ref_str, tx, pos, ref, alt);
            if (alt_str == NULL) return mc_exon_splice_sites;
            ret = predict_molecular_consequence(ref_str, alt_str, tx->is_mito);
            // Rprintf("alt : %s\n\n", alt_str);
            free(alt_str);
        }
    }
    return ret;
}
static struct dict *wnames = NULL;
struct csq *predict_func(const char *chr, int pos, int strand, const char *ref, const char *alt,
                         struct gtf_spec *G,
                         const char *fasta, faidx_t *fai, int *n, int debug)
{
    *n = 0;
    
    struct region_itr *itr = gtf_query(G, chr, pos - 1000, pos + 1000);
    if (itr == NULL) {
        int id = dict_query(G->name, chr);
        if (id == -1) {
            if (wnames == NULL) wnames = dict_init();
            int idx = dict_query(wnames, chr);
            if (idx == -1) {
                warnings("Chromosome %s not found in GTF, use wrong database? ", chr);
                dict_push(wnames, chr);
            }
        }
        return NULL;
    }

    int n_csq = 0;
    int m_csq = 4;
    struct csq *csq = malloc(4*sizeof(struct csq));

    int j;
    for (j = 0; j < itr->n; ++j) {
        // for each gene
        struct gtf *g = (struct gtf*)itr->rets[j];
        if (g->type != feature_gene) continue;

        if (pos < g->start) {
            if (pos > g->start - 1000) {
                if (n_csq == m_csq) {
                    m_csq = m_csq*2;
                    csq = realloc(csq, sizeof(struct csq)*m_csq);
                }
                struct csq *csq0 = &csq[n_csq++];
                if (g->strand == strand) {
                    csq0->csq = strand == BED_STRAND_FWD ? mc_upstream_1K : mc_downstream_1K;
                } else {
                    csq0->csq = strand == BED_STRAND_REV ? mc_antisense_upstream_1K : mc_antisense_downstream_1K;
                }
                csq0->gene_name = g->gene_name;
            }
        } else if (pos > g->end) {
            if (pos < g->end + 1000) {
                if (n_csq == m_csq) {
                    m_csq = m_csq*2;
                    csq = realloc(csq, sizeof(struct csq)*m_csq);
                }
                struct csq *csq0 = &csq[n_csq++];
                if (g->strand == strand) {
                    csq0->csq = strand == BED_STRAND_FWD ? mc_downstream_1K : mc_upstream_1K;
                } else {
                    csq0->csq = strand == BED_STRAND_REV ? mc_antisense_downstream_1K : mc_antisense_upstream_1K;
                }
                csq0->gene_name = g->gene_name;
            }
        } else {
            int k;
            for (k = 0; k < g->n_gtf; ++k) {
                struct gtf *tx = g->gtf[k];
                if (tx->type != feature_transcript) continue;

                if (pos < tx->start) continue;
                if (pos > tx->end) continue;
                
                if (n_csq == m_csq) {
                    m_csq = m_csq*2;
                    csq = realloc(csq, sizeof(struct csq)*m_csq);
                }
                struct csq *csq0 = &csq[n_csq++];
                csq0->csq = predict_func0(chr, pos, strand, ref, alt, G, tx, fasta, fai, debug);
                csq0->gene_name = g->gene_name;
                if (csq0->csq == mc_reference) break;
            }
        }
    }
    
    region_itr_destroy(itr);
    
    *n = n_csq;

    if (n_csq == 0) {
        free(csq);
        return NULL;
    }
    
    return csq;
}
SEXP anno_conseq(SEXP _chr, SEXP _pos, SEXP _ref, SEXP _alt, SEXP _strand, SEXP _db, SEXP _fa, SEXP _debug)
{
    int l = Rf_length(_chr);
    if (Rf_length(_pos) != l) {
        Rprintf("Inconsistance length of chr and start position.\n");
        return R_NilValue;
    }

    Rboolean debug = asLogical(_debug);

    struct gtf_spec *G = (struct gtf_spec *)R_ExternalPtrAddr(_db);
    const char *fasta = translateChar(STRING_ELT(_fa, 0));
    faidx_t *fai = fai_load(fasta);
    
    if (fai == NULL)
        return(mkString("Failed to load faidx of fasta, use `samtools faidx` index FASTA first."));

    memory_init();
    
    SEXP sl = PROTECT(allocVector(VECSXP, 2));
    SEXP gene = PROTECT(allocVector(STRSXP, l));
    SEXP conseq = PROTECT(allocVector(STRSXP, l));
    // struct dict *wnames = NULL;
               
    for (int i = 0; i < l; ++i) {
        const char *chr = translateChar(STRING_ELT(_chr, i));
        int pos = INTEGER(_pos)[i];
        const char *ref = translateChar(STRING_ELT(_ref, i));
        const char *alt = translateChar(STRING_ELT(_alt, i));
        
        int strand = -1;
        if (!Rf_isNull(_strand)) {
            const char *s = translateChar(STRING_ELT(_strand, i));
            if (s[0] == '+') strand = BED_STRAND_FWD;
            else if (s[0] == '-') strand = BED_STRAND_REV;
        }

        int n_csq = 0;
        struct csq *csq = predict_func(chr, pos, strand, ref, alt, G, fasta, fai, &n_csq, debug);
        if (n_csq > 1) qsort(csq, n_csq, sizeof(struct csq), cmp);
        if (n_csq == 0) {
            SET_STRING_ELT(gene, i, mkChar("."));
            SET_STRING_ELT(conseq, i, mkChar("intergenic"));
        } else {
            SET_STRING_ELT(gene, i, mkChar(GTF_genename(G, csq[0].gene_name)));
            SET_STRING_ELT(conseq, i, mkChar(MCT[csq[0].csq]));
            free(csq);
        }
    }

    dict_destroy(wnames);
    fai_destroy(fai);

    memory_release();
    
    SET_VECTOR_ELT(sl, 0, gene);
    SET_VECTOR_ELT(sl, 1, conseq);
    UNPROTECT(3);
    return sl;
}


