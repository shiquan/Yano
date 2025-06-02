#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Matrix.h>
#include "utils.h"
#include "dict.h"
#include "perm.h"
#include "variant_type.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"
#include "htslib/bgzf.h"
#include "htslib/faidx.h"

SEXP vcf_sample_names(SEXP _vcf)
{
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

    bcf_hdr_t *hdr = bcf_hdr_read(fp);    
    if (hdr == NULL) {
        Rprintf("Failed to load header of %s.\n", vcf_fname);
        return R_NilValue;
    }

    int nsample = bcf_hdr_nsamples(hdr);
    SEXP sm = PROTECT(allocVector(STRSXP, nsample));
    int i;
    for (i = 0; i < nsample; ++i) {
        const char *smp = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
        SET_STRING_ELT(sm, i, mkChar(smp));
    }
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    UNPROTECT(1);
    return sm;
}
SEXP vcf_query_value(SEXP _vcf, SEXP _chr, SEXP _pos, SEXP _ref)
{
    const char *vcf_fname = translateChar(STRING_ELT(_vcf, 0));
    const char *chr = translateChar(STRING_ELT(_chr, 0));
    char *ref = NULL;
    if (_ref != R_NilValue) {
        ref = (char*)translateChar(STRING_ELT(_ref, 0));
    }
    
    int pos = asInteger(_pos);   
    htsFile *fp = hts_open(vcf_fname, "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    int rid = bcf_hdr_name2id(hdr, chr);
    if (rid == -1) {
        hts_close(fp);
        bcf_hdr_destroy(hdr);
        error("No such chromosome, %s", chr);
    }
    hts_idx_t *idx = bcf_index_load(vcf_fname);
    if (idx == NULL) {
        hts_close(fp);
        bcf_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        error("No index file, plz index the VCF/BCF file first.");
    }

    hts_itr_t *itr = bcf_itr_queryi(idx, rid, pos-1, pos);
    if (!itr) {
        hts_close(fp);
        bcf_hdr_destroy(hdr);
        hts_idx_destroy(idx);
        error("No such position, %s, %d", chr, pos);        
    }
    bcf1_t *v = bcf_init();
    int ret;
    for (;;) {
        ret = bcf_itr_next(fp, itr, v);
        if (ret < 0) break;
        if (v->pos +1 == pos) {
            if (ref) {
                bcf_unpack(v, BCF_UN_INFO |BCF_UN_FMT);
                if (strcmp(ref, v->d.allele[0]) == 0) break;
                else continue;
            } else {
                break;
            }
        }
        if (v->pos > pos) {
            ret = -1;
            break;   
        }
    }
    if (ret < 0) {
        bcf_destroy(v);
        hts_itr_destroy(itr);
        hts_idx_destroy(idx);
        hts_close(fp);
        bcf_hdr_destroy(hdr);
        return R_NilValue;
    }

    bcf_unpack(v, BCF_UN_INFO |BCF_UN_FMT);

    int32_t *gt_arr = NULL;
    int ngt_arr = 0;
    int ngt = bcf_get_genotypes(hdr, v, &gt_arr, &ngt_arr);
    int nsample = bcf_hdr_nsamples(hdr);
    int max_ploidy = ngt/nsample;

    SEXP df = PROTECT(allocVector(VECSXP, max_ploidy));
    SEXP ta[max_ploidy];
    
    int i, j;
    for (i = 0; i < max_ploidy; ++i) {
        ta[i] = PROTECT(allocVector(INTSXP, nsample));
        for (j = 0; j < nsample; ++j) {
            INTEGER(ta[i])[j] = 0;
        }
    }
    
    for (i = 0; i < nsample; i++) {
        int32_t *ptr = gt_arr + i*max_ploidy;
        for (j = 0; j < max_ploidy; ++j) {
            if (ptr[j] == bcf_int32_vector_end) break;
            if (bcf_gt_is_missing(ptr[j])) continue;
            int allele = bcf_gt_allele(ptr[j]);
            INTEGER(ta[allele])[i] = INTEGER(ta[allele])[i] +1 ;
        }
    }

    if (!ref) {
        // warnings        
    }    

    for (i = 0; i < max_ploidy; ++i) {
        SET_VECTOR_ELT(df, i, ta[i]);
    }

    setAttrib(df, R_ClassSymbol, mkString("data.frame"));
    SEXP names = PROTECT(allocVector(STRSXP, max_ploidy));
    for (i = 0; i < max_ploidy; ++i) {
        SET_STRING_ELT(names,  i, mkChar(v->d.allele[i]));
    }

    setAttrib(df, R_NamesSymbol, names);

    SEXP rownames = PROTECT(allocVector(STRSXP, nsample));
    for (i = 0; i < nsample; ++i) {
        const char *smp = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
        SET_STRING_ELT(rownames,  i, mkChar(smp));
    }

    setAttrib(df, R_RowNamesSymbol, rownames);

    free(gt_arr);
    bcf_hdr_destroy(hdr);
    bcf_destroy(v);
    hts_itr_destroy(itr);
    hts_idx_destroy(idx);
    hts_close(fp);
    UNPROTECT(max_ploidy+3);
    return df;
}

extern void smooth_W(double * const a, double *s, const int N, CHM_SP W);
extern void smooth_W2(double * const a, double *s, const int N, CHM_SP W);

struct ret {
    int rid;
    int pos;
    char *ref;
    char *alt;
    double ts; // t value
};

static int n_ret = 0;
static int m_ret = 0;
static struct ret **rets = NULL;

void buffer_init() {
    m_ret = 1000;
    rets = malloc(m_ret*sizeof(struct ret*));
}
static void write_out(struct ret *ret)
{
    if (n_ret == m_ret) {
        m_ret = m_ret *2;
        rets = realloc(rets, m_ret*sizeof(struct ret*));
    }
    if (ret) rets[n_ret++] = ret;
}
struct ret *ret_create() {
    struct ret *ret = malloc(sizeof(struct ret));
    memset(ret, 0, sizeof(struct ret));
    return ret;
}

static bcf1_t *read_line(htsFile *fp, bcf_hdr_t *hdr)
{
    bcf1_t *v = bcf_init();
    for (;;) {
        int ret = bcf_read(fp, hdr, v);
        if (ret != 0) {
            bcf_destroy(v);
            v = NULL;
            break;
        }
        if (v->rid == -1) continue;
        break;
    }
    return v;
}

static struct ret *run_it(bcf_hdr_t *hdr, bcf1_t *v, int n_unit, struct dict *samples, int perm, double maf, CHM_SP W)
{
    bcf_unpack(v, BCF_UN_INFO);

    int32_t *gt_arr = NULL;
    int ngt_arr = 0;
    int ngt = bcf_get_genotypes(hdr, v, &gt_arr, &ngt_arr);
    if (ngt < 1) {
        free(gt_arr);
        return NULL;
    }
    int nsample = bcf_hdr_nsamples(hdr);
    int max_ploidy = ngt/nsample;
    if (max_ploidy < 2) {
        free(gt_arr);
        return NULL;
    }

    int ac[2]; // ref and nonref
    double a0[n_unit]; // ref count per sample
    double a1[n_unit];
    ac[0] = 0; ac[1] = 0;
    // filtering
    memset(a0, 0, n_unit*sizeof(double));
    memset(a1, 0, n_unit*sizeof(double));
    int i, j;        
    for (i = 0; i < nsample; i++) {
        int idx = dict_query(samples, hdr->samples[i]);
        if (idx < 0) continue;
        int32_t *ptr = gt_arr + i*max_ploidy;
        for (j = 0; j < max_ploidy; ++j) {
            if (ptr[j] == bcf_int32_vector_end) break;
            if (bcf_gt_is_missing(ptr[j])) continue;
            
            int allele = bcf_gt_allele(ptr[j]);
            if (allele) {
                ac[1]++;
                a1[idx]++;
            } else {
                ac[0]++;
                a0[idx]++;
            }
        }
    }
    if ((double)ac[1]/(ac[0]+ac[1]) < maf || (double)ac[0]/(ac[0]+ac[1]) < maf) {
        free(gt_arr);
        return NULL;
    }

    struct ret *ret = ret_create();
    ret->rid = v->rid;
    ret->pos = v->pos+1;
    
    double tmpa[n_unit];
    double tmpb[n_unit];
    double tmpa_s[n_unit];
    double tmpb_s[n_unit];
    
    // only test for nonref alleles
    // SPA

    smooth_W2(a1, tmpa, n_unit, W);
    smooth_W2(a0, tmpb, n_unit, W);

    /* for (j = 0; j < n_unit; ++j) {         */
    /*     tmpa[j] = a1[j]; */
    /*     tmpb[j] = a0[j]; */
    /* }             */

    smooth_W(tmpa, tmpa_s, n_unit, W);
    smooth_W(tmpb, tmpb_s, n_unit, W);
    double mna = 0, mnb = 0;
    double mna_s = 0, mnb_s = 0;
    for (j = 0; j < n_unit; ++j) {
        mna += tmpa[j];
        mnb += tmpb[j];
        mna_s += tmpa_s[j];
        mnb_s += tmpb_s[j];
    }
    mna_s = mna_s/n_unit;
    mnb_s = mnb_s/n_unit;
    mna = mna/n_unit;
    mnb = mnb/n_unit;
    
    double Lx1 = 0, Lx2 = 0; // Ly1 = 0, Ly2 = 0;
    double ra = 0, rb1 = 0, rb2 = 0;
    for (j = 0; j < n_unit; ++j) {
        Lx1 += pow(tmpa_s[j]-mna,2);
        Lx2 += pow(tmpa[j]-mna,2);
        //Ly1 += pow(tmpb_s[j]-mnb,2);
        //Ly2 += pow(tmpb[j]-mnb,2);
        tmpa_s[j] = tmpa_s[j] - mna_s;
        tmpb_s[j] = tmpb_s[j] - mnb_s;
        ra  += tmpa_s[j] * tmpb_s[j];
        rb1 += pow(tmpa_s[j],2);
        rb2 += pow(tmpb_s[j],2);
    }
    
    rb1 = sqrt(rb1);
    rb2 = sqrt(rb2);
    
    double Lx = Lx1/Lx2;
    //double Ly = Ly1/Ly2;
    double r = ra/(rb1*rb2);
    double e = sqrt(Lx) * (1-r);
    double es[perm];
    int k;
    double mean = 0;
    for (k = 0; k < perm; ++k) {
        shuffle_index(tmpa, k, n_unit);
        smooth_W(tmpa, tmpa_s, n_unit, W);
        mna_s = 0;
        for (j = 0; j < n_unit; ++j) {
            mna_s += tmpa_s[j];
        }
        mna_s = mna_s/n_unit;
        Lx1 = 0;
        ra = 0;
        rb1 = 0;
        for (j = 0; j < n_unit; ++j) {
            Lx1 += pow(tmpa_s[j]-mna, 2);
            tmpa_s[j] = tmpa_s[j] - mna_s;
            rb1 += pow(tmpa_s[j], 2);
        }
        rb1 = sqrt(rb1);
        Lx = Lx1/Lx2;
        
        for (j = 0; j < n_unit; ++j) {
            ra += tmpa_s[j] * tmpb_s[j];
        }
        
        r = ra/(rb1*rb2);
        es[k] = sqrt(Lx) *(1-r);
        mean += es[k];
    }
    mean = mean/perm;
    double var = 0;
    for (j = 0; j < perm; ++j) {
        var += pow((es[j] -mean),2);
    }
    var = sqrt(var/perm);
    ret->ts = (e - mean)/var;
    
    ret->ref = strdup(v->d.allele[0]);
    if (max_ploidy == 2) {
        ret->alt = strdup(v->d.allele[1]);
    } else {
        kstring_t tmp = {0,0,0};
        int k;
        for (k = 1; k < max_ploidy; ++k) {
            if (k > 1) kputc(',', &tmp);
            kputs(v->d.allele[k], &tmp);
        }
        ret->alt = tmp.s;
    }
    
    free(gt_arr);

    return ret;
}

SEXP spa_vcf(SEXP _vcf, SEXP _W, SEXP _maf, SEXP _permut, SEXP _threads)
{
    const int perm = asInteger(_permut);
    int n_thread = asInteger(_threads);

    CHM_SP W = AS_CHM_SP__(_W);
    const char *vcf_fname = translateChar(STRING_ELT(_vcf, 0));
    double maf = asReal(_maf);    
    htsFile *fp = hts_open(vcf_fname, "r");
    if (fp == NULL) {
        Rprintf("Failed to read %s.\n", vcf_fname);
        return R_NilValue;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);    

    struct dict *samples = dict_init();
    int n_unit = W->nrow;
    SEXP Matrix_DimNamesSym = install("Dimnames");
    SEXP dn = PROTECT(GET_SLOT(_W, Matrix_DimNamesSym));
    SEXP rn = VECTOR_ELT(dn, 0);
    if (isNull(rn)) error("No row names.");
    UNPROTECT(1);
    
    int i;
    for (i = 0; i < n_unit; ++i) {
        dict_push(samples, CHAR(STRING_ELT(rn, i)));
    }

    int nsample = bcf_hdr_nsamples(hdr);

    if (nsample < W->nrow)
        error("Inconsistant sample names in the weight matrix and the VCF header.");

    random_index_init(perm, n_unit);
    
#pragma omp parallel num_threads(n_thread)
    for (;;) {
        bcf1_t *v = NULL;
#pragma omp critical (read)
        v = read_line(fp, hdr);
        if (v == NULL) break;

        struct ret *ret = run_it(hdr, v, n_unit, samples, perm, maf, W);
        if (ret == NULL) continue;
        
#pragma omp critical (write)
        write_out(ret);        
    }

    if (n_ret == 0) {
        dict_destroy(samples);
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        random_index_free();
        free(rets);
        return R_NilValue;
    }
    
    SEXP Vchr = PROTECT(allocVector(STRSXP, n_ret));
    SEXP Vpos = PROTECT(allocVector(INTSXP, n_ret));
    SEXP Vref = PROTECT(allocVector(STRSXP, n_ret));
    SEXP Valt = PROTECT(allocVector(STRSXP, n_ret));
    SEXP Vt = PROTECT(allocVector(REALSXP, n_ret));

    for (i = 0; i < n_ret; ++i) {
        struct ret *ret = rets[i];
        const char *chr = hdr->id[BCF_DT_CTG][ret->rid].key;
        SET_STRING_ELT(Vchr, i, mkChar(chr));
        INTEGER(Vpos)[i] = ret->pos;
        if (ret->ref) {
            SET_STRING_ELT(Vref, i, mkChar(ret->ref));
        }

        if (ret->alt) {
            SET_STRING_ELT(Valt, i, mkChar(ret->alt));
        }
        
        REAL(Vt)[i] = ret->ts;
    }
    
    dict_destroy(samples);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    random_index_free();

    for (i = 0; i < n_ret; ++i) {
        struct ret *ret = rets[i];
        if (ret->ref) free(ret->ref);
        if (ret->alt) free(ret->alt);
        free(ret);
    }
    free(rets);

    SEXP ta = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ta, 0, Vchr);
    SET_VECTOR_ELT(ta, 1, Vpos);
    SET_VECTOR_ELT(ta, 2, Vref);
    SET_VECTOR_ELT(ta, 3, Valt);
    SET_VECTOR_ELT(ta, 4, Vt);

    setAttrib(ta, R_ClassSymbol, mkString("data.frame"));
    SEXP names = PROTECT(allocVector(STRSXP, 5));
    SET_STRING_ELT(names,  0, mkChar("chr"));
    SET_STRING_ELT(names,  1, mkChar("pos"));
    SET_STRING_ELT(names,  2, mkChar("ref"));
    SET_STRING_ELT(names,  3, mkChar("alt"));
    SET_STRING_ELT(names,  4, mkChar("tval"));    
    setAttrib(ta, R_NamesSymbol, names);

    SEXP rownames = PROTECT(allocVector(INTSXP, 2));
    SET_INTEGER_ELT(rownames, 0, NA_INTEGER);
    SET_INTEGER_ELT(rownames, 1, -n_ret);
    setAttrib(ta, R_RowNamesSymbol, rownames);
  
    UNPROTECT(8);

    return ta;
}
