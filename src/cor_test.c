#include <omp.h>

#include<R.h>
#include<Rdefines.h>
//#include<Rinternals.h>

#ifndef Matrix_stubs
#define Matrix_stubs 1
#include "Matrix_stubs.c"
#endif

#include <assert.h>


void progress(int p) {
    static int displayed = -1;
    static int count = 0;
    static char bar[] = "================================================== ";

    count++;

    int r = count*50;
    
    if (displayed==-1) {
        if (r>50) return;
        REprintf("|--------------------------------------------------|\n|");            
        displayed = 0;
    }

    int toPrint = r-displayed;
    if (toPrint==0) return;
    bar[toPrint] = '\0';
    REprintf("%s", bar);
    bar[toPrint] = '=';
    displayed = p;
    if (p==50) {
        REprintf("|\n");
        displayed = -1;
    }
    R_FlushConsole();
}
// # nocov end
static cholmod_common c;

#define CACHE_PER_BATCH 10000
#define MIN_HIT 1

void shuffle_double_arr(double *a, const int n)
{
    int i;
    for (i = 0; i < n-1; ++i) {
        int j = i + rand() / (RAND_MAX / (n - i) + 1);
        double t = a[j];
        a[j] = a[i];
        a[i] = t;
    }
}
int *random_idx(const int n)
{
    int *idx = R_Calloc(n, int);
    int i;
    for (i = 0; i < n; ++i) idx[i] = i;
    for (i = 0; i < n-1; ++i) {
        int j = i + rand() / (RAND_MAX / (n - i) + 1);
        int t = idx[j];
        idx[j] = idx[i];
        idx[i] = t;
    }
    return idx;
}
static void shuffle(double tmp[], int const idx[], const int n)
{
    double aux[n];
    int i;
    for (i = 0; i < n; i++) aux[idx[i]] = tmp[i];
    for (i = 0; i < n; i++) tmp[i] = aux[i];
}

SEXP cor_test(SEXP _A, SEXP _B, SEXP _trans)
{
    Rboolean tr = asLogical(_trans);
    
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP B = AS_CHM_SP__(_B);
    R_CheckStack();

    if (A->ncol != B->ncol || A->nrow != B->nrow) return mkString("Unequal matrix size.");
    if (A->ncol < 2) return mkString("Too few cells.");

    if (tr) {
        A = M_cholmod_transpose(A, (int)A->xtype, &c);
        B = M_cholmod_transpose(B, (int)B->xtype, &c);
    }
    
    int N = A->nrow;
    int nc = A->ncol;

    SEXP rval = PROTECT(allocVector(REALSXP, nc));
    
    double *tmpa = R_Calloc(N, double);
    double *tmpb = R_Calloc(N, double);

    int *ap = (int*)A->p;
    int *ai = (int*)A->i;
    double *ax = (double*)A->x;

    int *bp = (int*)B->p;
    int *bi = (int*)B->i;
    double *bx = (double*)B->x;

    int ci;
    for (ci = 0; ci < nc; ++ci) {
        if (ap[ci] == ap[ci+1]) {
            REAL(rval)[ci] = 0;
            continue;   
        }
        if (bp[ci] == bp[ci+1]) {
            REAL(rval)[ci] = 0;
            continue;               
        }
        
        memset(tmpa, 0, sizeof(double)*N);
        memset(tmpb, 0, sizeof(double)*N);
        // scale
        int j;
        double mn = 0.0;
        //double sd = 0.0;

        for (j = ap[ci]; j < ap[ci+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            mn += ax[j];
            tmpa[ai[j]] = ax[j];
        }
        mn = mn/N;

        for (j = 0; j < N; ++j) tmpa[j] = tmpa[j] - mn;
        
        mn = 0;
        //sd = 0;        
        for (j = bp[ci]; j < bp[ci+1]; ++j) {
            if (ISNAN(bx[j])) continue;
            mn += bx[j];
            tmpb[bi[j]] = bx[j];
        }
        mn = mn/N;

        for (j = 0; j < N; ++j) tmpb[j] = tmpb[j] - mn;

        double r = 0,
            ra = 0,
            rb = 0;
        
        for (j = 0; j < N; ++j) {
            r += tmpa[j]*tmpb[j];
            ra += pow(tmpa[j],2);
            rb += pow(tmpb[j],2);
        }

        REAL(rval)[ci] = r/ (sqrt(ra) * sqrt(rb));
    }

    R_Free(tmpa);
    R_Free(tmpb);

    if (tr) {
        M_cholmod_free_sparse(&A, &c);
        M_cholmod_free_sparse(&B, &c);
    }

    UNPROTECT(1);
    return rval;
}

void smooth_W(double * const a, double *s, const int N, CHM_SP W)
{
    int *wp = (int*)W->p;
    int *wi = (int*)W->i;
    double *wx = (double*)W->x;

    memset(s, 0, N*sizeof(double));
    
    int i;
    for (i = 0; i < N; ++i) {
        if (wp[i] == wp[i+1]) continue;
        int j;
        for (j = wp[i]; j < wp[i+1]; ++j) {
            int idx = wi[j];
            s[i] += wx[j] * a[idx];
        }
    }
}

SEXP moransi_mc_test(SEXP _A, SEXP _W, SEXP _trans, SEXP _permut, SEXP _threads)
{
    Rboolean tr = asLogical(_trans);
    int perm = asInteger(_permut);
    int n_thread = asInteger(_threads);
    
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP W = AS_CHM_SP__(_W);
    
    if (A->stype) return mkString("A cannot be symmetric");
    if (W->stype) return mkString("W cannot be symmetric");
    
    //double one[] = {1, 0};

    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (A->ncol < 2) return mkString("Too few cells."); // to do

    if (tr) {
        A = M_cholmod_transpose(A, (int)A->xtype, &c);
    }

    R_CheckStack();
    
    int N_cell = A->nrow;
    int N_feature = A->ncol;

    SEXP Ival = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Pval = PROTECT(allocVector(REALSXP, N_feature));    

    const int *ap = (int*)A->p;
    const int *ai = (int*)A->i;
    const double *ax = (double*)A->x;

    const int *wp = (int*)W->p;
    const int *wi = (int*)W->i;
    const double *wx = (double*)W->x;

    int *hits = R_Calloc(N_feature, int);    
    int *steps = R_Calloc(N_feature, int);
    int N = N_feature;
    for (int i = 0; i < N; ++i) steps[i] = i;    
    int n_step = perm/CACHE_PER_BATCH;
    if (n_step < 1) n_step = 1;
    
    for (int step = 0; step < n_step; ++step) {
        int **ris = NULL;
        if (perm > 1) {
            ris = R_Calloc(CACHE_PER_BATCH,int*);
            for (int pi = 0; pi < CACHE_PER_BATCH; ++pi) {
                ris[pi] = random_idx(N_cell);
            }
        }

        int i;
#pragma omp parallel for num_threads(n_thread)
        for (i = 0; i < N; ++i) {
            int idx = steps[i];
            if (ap[idx] == ap[idx+1]) {
                REAL(Ival)[idx] = 0;
                REAL(Pval)[idx] = 0;
                continue;
            }
            
            double *tmp = R_Calloc(N_cell, double);
            int j;
            double mn = 0,
                //sd = 0,
                xix = 0,
                xij = 0;
            
            for (j = ap[idx]; j < ap[idx+1]; ++j) {
                if (ISNAN(ax[j])) continue;
                mn += ax[j];
                tmp[ai[j]] = ax[j];
            }
            mn = mn/N_cell;
            
            for (j = 0; j < N_cell; ++j) {
                tmp[j] = tmp[j] - mn;
                xix += pow(tmp[j], 2);
            }
            
            for (int k = 0; k < N_cell; ++k) {
                for (j = wp[k]; j < wp[k+1]; ++j) {
                    if (ISNAN(wx[j])) continue;
                    int cid = wi[j];
                    if (k == cid) continue; // i != j
                    xij += wx[j] * tmp[k] * tmp[cid];
                }
            }
            double I = xij/xix;
            // double p = 0;
            
            if (perm > 1) {
                double *Is = R_Calloc(CACHE_PER_BATCH, double);
                xij = 0;
                
                for (int pi = 0; pi < CACHE_PER_BATCH; ++pi) {
                    shuffle(tmp, ris[pi], N_cell);
                    for (int k = 0; k < N_cell; ++k) {
                        for (j = wp[k]; j < wp[k+1]; ++j) {
                            if (ISNAN(wx[j])) continue;
                            int cid = wi[j];
                            if (k == cid) continue; // i != j
                            xij += wx[j] * tmp[k] * tmp[cid];
                        }
                    }
                    Is[pi] = xij/xix;
                }
                for (int k = 0; k < CACHE_PER_BATCH; ++k) {
                    //Rprintf("%f\t%f\n",e, es[k]);
                    if (Is[k]>=I) hits[idx] += 1;
                }
                // p = p/(double)perm;
                R_Free(Is);
            }
            R_Free(tmp);

#pragma omp critical
            if (step == 1) {            
                REAL(Ival)[i] = I;
                // REAL(Pval)[i] = p;
            }
        }

        if (perm > 1) {
            for (int pi = 0; pi < CACHE_PER_BATCH; ++pi) R_Free(ris[pi]);
            R_Free(ris);
        }

        int j = 0;
        for (i = 0; i < N_feature; ++i) {
            if (hits[i] == -1) continue;
            if (hits[i] > MIN_HIT) {
                REAL(Pval)[i] = (double)hits[i]/((step+1)*CACHE_PER_BATCH);
            } else {
                steps[j++] = i;
            }
        }
        N = j;
    }

    for (int i = 0; i < N; ++i) {
        int idx = steps[i];
        REAL(Pval)[idx] = (double)hits[idx]/(n_step*CACHE_PER_BATCH);
    }
    
    R_Free(hits);
    R_Free(steps);
    
    if (tr) {
        M_cholmod_free_sparse(&A, &c);
    }

    SEXP ta = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ta, 0, Ival);
    SET_VECTOR_ELT(ta, 1, Pval);

    UNPROTECT(3);
    return ta;
}
 
SEXP smooth_test(SEXP _A, SEXP _W,
                 SEXP idx,
                 SEXP cs, SEXP _scale)
            
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP W = AS_CHM_SP__(_W);

    const int scale_factor = asInteger(_scale);

    if (A->stype) return mkString("A cannot be symmetric");
    if (W->stype) return mkString("W cannot be symmetric");
        
    //double one[] = {1, 0};

    // if (A->ncol != B->ncol || A->nrow != B->nrow) return mkString("A and B do not match");
    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (A->ncol < 2) return mkString("Too few cells."); // to do

    A = M_cholmod_transpose(A, (int)A->xtype, &c);
    W = M_cholmod_transpose(W, (int)W->xtype, &c);

    R_CheckStack();
    const int N_cell = A->nrow;
    const int N_feature = length(idx);

    const int * ap = (int*)A->p;
    const int * ai = (int*)A->i;
    const double * ax = (double*)A->x;

    SEXP ta = PROTECT(allocVector(VECSXP, N_feature));

    int i;
    for (i = 0; i < N_feature; ++i) {
        int ii = INTEGER(idx)[i]  -1;
        if (ap[ii] == ap[ii+1]) {
            continue;
        }
        
        double *tmpa = R_Calloc(N_cell, double);
        memset(tmpa, 0, sizeof(double)*N_cell);
        int j;
        for (j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            int cid = ai[j];
            tmpa[cid] = log(ax[j]/REAL(cs)[cid]*scale_factor + 1);
        }
        
        double *tmpa_s = R_Calloc(N_cell, double);
        smooth_W(tmpa, tmpa_s, N_cell, W);
        SEXP val = PROTECT(allocVector(REALSXP, N_cell));
        for (j = 0; j < N_cell; ++j) {
            REAL(val)[j] = tmpa_s[j];
        }
        SET_VECTOR_ELT(ta, i, val);
        
        R_Free(tmpa_s);
        R_Free(tmpa);
    }

    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&W, &c);

    UNPROTECT(N_feature+1);
    return ta;
}
SEXP D_test(SEXP _A, SEXP _B, SEXP _W,
            SEXP _permut,
            SEXP _threads,
            SEXP idx, SEXP bidx,
            //SEXP cidx,
            SEXP cs, SEXP _scale,
            SEXP _sens)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP B = AS_CHM_SP__(_B);
    CHM_SP W = AS_CHM_SP__(_W);

    const int perm = asInteger(_permut);
    const int n_thread = asInteger(_threads);
    const int scale_factor = asInteger(_scale);

    Rboolean sens = asLogical(_sens);
    
    if (A->stype) return mkString("A cannot be symmetric");
    if (B->stype) return mkString("B cannot be symmetric");
    if (W->stype) return mkString("W cannot be symmetric");
    
    // double one[] = {1, 0};

    // if (A->ncol != B->ncol || A->nrow != B->nrow) return mkString("A and B do not match");
    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (A->ncol < 2) return mkString("Too few cells."); // to do

    A = M_cholmod_transpose(A, (int)A->xtype, &c);
    B = M_cholmod_transpose(B, (int)B->xtype, &c);
    W = M_cholmod_transpose(W, (int)W->xtype, &c);

    R_CheckStack();
    const int N_cell = A->nrow;
    //const int N_cell = length(cidx);
    const int N_feature = length(idx);

    assert (length(bidx) == N_feature);
    
    SEXP LXval = PROTECT(allocVector(REALSXP, N_feature));
    SEXP LYval = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Rval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Dval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Tval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Mval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Vval  = PROTECT(allocVector(REALSXP, N_feature));
    
    const int * const ap = (int*)A->p;
    const int * const ai = (int*)A->i;
    const double * const ax = (double*)A->x;

    const int *const bp = (int*)B->p;
    const int *const bi = (int*)B->i;
    const double *const bx = (double*)B->x;
    
    int **ris = NULL;
    ris = R_Calloc(perm,int*);
    for (int pi = 0; pi < perm; ++pi) {
        ris[pi] = random_idx(N_cell);
    }
    
    int i;
    
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < N_feature; ++i) {
        // todo: bidx
        int ii = INTEGER(idx)[i]  -1;
        int ij = INTEGER(bidx)[i] -1;
        //Rprintf("%d\t%d\n", ii, ij);
        if (ap[ii] == ap[ii+1] || bp[ij] == bp[ij+1]) {
            REAL(LXval)[i] = 0;
            REAL(LYval)[i] = 0;
            REAL(Rval)[i]  = 0;
            REAL(Dval)[i]  = 0;
            REAL(Tval)[i]  = 0;
            continue;
        }
        
        double *tmpa = R_Calloc(N_cell, double);
        double *tmpb = R_Calloc(N_cell, double);
        memset(tmpa, 0, sizeof(double)*N_cell);
        memset(tmpb, 0, sizeof(double)*N_cell);
        double mna = 0,
            mnb = 0;
        int j;
        for (j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            int cid = ai[j];
            //cid = INTEGER(cidx)[cid]-1;
            tmpa[cid] = ax[j];
        }
        
        for (j = bp[ij]; j < bp[ij+1]; ++j) {
            if (ISNAN(bx[j])) continue;
            int cid = bi[j];
            //cid = INTEGER(cidx)[cid]-1;
            tmpb[cid] = bx[j];
            if (sens) {
                tmpb[cid] = tmpb[cid] - tmpa[cid];
                if (tmpb[cid] < 0) tmpb[cid] = 0;
            }

            tmpb[cid] = log(tmpb[cid]/REAL(cs)[cid]*scale_factor + 1);
            mnb += tmpb[cid];
        }
        mnb = mnb/N_cell;

        for (j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;            
            int cid = ai[j];
            //cid = INTEGER(cidx)[cid]-1;
            tmpa[cid] = log(ax[j]/REAL(cs)[cid]*scale_factor + 1);
            mna += tmpa[cid];
        }
        mna = mna/N_cell;
        
        double *tmpa_s = R_Calloc(N_cell, double);
        double *tmpb_s = R_Calloc(N_cell, double);
        smooth_W(tmpa, tmpa_s, N_cell, W);
        smooth_W(tmpb, tmpb_s, N_cell, W);
        double mna_s = 0,
            mnb_s = 0;
        for (j = 0; j < N_cell; ++j) {
            mna_s += tmpa_s[j];
            mnb_s += tmpb_s[j];
        }
        mna_s = mna_s/(double)N_cell;
        mnb_s = mnb_s/(double)N_cell;
        
        double Lx1 = 0,
            Lx2 = 0,
            Ly1 = 0,
            Ly2 = 0,
            ra  = 0,
            rb1 = 0,
            rb2 = 0;
        for (j = 0; j < N_cell; ++j) {
            Lx1 += pow(tmpa_s[j]-mna,2);
            Lx2 += pow(tmpa[j]-mna,2);
            Ly1 += pow(tmpb_s[j]-mnb,2);
            Ly2 += pow(tmpb[j]-mnb,2);

            tmpa_s[j] = tmpa_s[j] - mna_s;
            tmpb_s[j] = tmpb_s[j] - mnb_s;
            ra  += tmpa_s[j] * tmpb_s[j];
            rb1 += pow(tmpa_s[j],2);
            rb2 += pow(tmpb_s[j],2);
        }
        
        rb1 = sqrt(rb1);
        rb2 = sqrt(rb2);
        
        double Lx = Lx1/Lx2;
        double Ly = Ly1/Ly2;
        double r = ra/(rb1*rb2);
        double e = sqrt(Lx) * (1-r);
#pragma omp critical
        {
            REAL(LXval)[i] = Lx;
            REAL(LYval)[i] = Ly;
            REAL(Rval)[i]  = r;
            REAL(Dval)[i]  = e;
        }

        // mean, var
        double mean = 0,
            var = 0;
        
        double *es = R_Calloc(perm, double);
        int k;
        for (k = 0; k < perm; ++k) {
            shuffle(tmpa, ris[k], N_cell);
            smooth_W(tmpa, tmpa_s, N_cell, W);
            mna_s = 0;
            for (j = 0; j < N_cell; ++j) {
                mna_s += tmpa_s[j];
            }
            mna_s = mna_s/(double)N_cell;
            Lx1 = 0;
            ra = 0;
            rb1 = 0;
            for (j = 0; j < N_cell; ++j) {
                Lx1 += pow(tmpa_s[j]-mna,2);
                tmpa_s[j] = tmpa_s[j] - mna_s;
                ra += tmpa_s[j] * tmpb_s[j];
                rb1 += pow(tmpa_s[j],2);
            }
            rb1 = sqrt(rb1);
            Lx = Lx1/Lx2;
            r = ra/(rb1*rb2);
            es[k] = sqrt(Lx) *(1-r);
            mean += es[k];
        }
        mean = mean/perm;

        for (k = 0; k < perm; ++k) {
            var += pow((es[k] -mean),2);
        }
        var = sqrt(var/perm);

        double t = (e - mean)/var;
        
        R_Free(es);
        R_Free(tmpa_s);
        R_Free(tmpb_s);
        R_Free(tmpa);
        R_Free(tmpb);
        
#pragma omp critical
        {
            REAL(Tval)[i]  = t;
            REAL(Mval)[i]  = mean;
            REAL(Vval)[i]  = var;
        }
    }
    
    for (int pi = 0; pi < perm; ++pi) R_Free(ris[pi]);
    R_Free(ris);

    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&B, &c);
    M_cholmod_free_sparse(&W, &c);

    SEXP ta = PROTECT(allocVector(VECSXP, 7));
    SET_VECTOR_ELT(ta, 0, LXval);
    SET_VECTOR_ELT(ta, 1, LYval);
    SET_VECTOR_ELT(ta, 2, Rval);
    SET_VECTOR_ELT(ta, 3, Dval);
    SET_VECTOR_ELT(ta, 4, Tval);
    SET_VECTOR_ELT(ta, 5, Mval);
    SET_VECTOR_ELT(ta, 6, Vval);

    UNPROTECT(8);
    return ta;
}
SEXP D_test_cell(SEXP _A, SEXP _B, SEXP _W,
                 SEXP _permut,
                 SEXP _threads,
                 SEXP idx, 
                 SEXP cs, SEXP _scale,
                 SEXP _sens)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP W = AS_CHM_SP__(_W);

    const int perm = asInteger(_permut);
    const int n_thread = asInteger(_threads);
    const int scale_factor = asInteger(_scale);

    Rboolean sens = asLogical(_sens);
    
    if (A->stype) return mkString("A cannot be symmetric");
    if (W->stype) return mkString("W cannot be symmetric");
    
    //double one[] = {1, 0};
    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (A->ncol < 2) return mkString("Too few cells."); // to do

    A = M_cholmod_transpose(A, (int)A->xtype, &c);
    W = M_cholmod_transpose(W, (int)W->xtype, &c);

    R_CheckStack();
    const int N_cell = A->nrow;
    const int N_feature = length(idx);
    
    SEXP LXval = PROTECT(allocVector(REALSXP, N_feature));
    // SEXP LYval = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Rval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Dval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Tval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Mval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Vval  = PROTECT(allocVector(REALSXP, N_feature));
    
    const int * ap = (int*)A->p;
    const int * ai = (int*)A->i;
    const double * ax = (double*)A->x;
    
    int **ris = NULL;
    ris = R_Calloc(perm,int*);
    for (int pi = 0; pi < perm; ++pi) {
        ris[pi] = random_idx(N_cell);
    }
    
    double *bx = malloc(N_cell*sizeof(double));
    for (int ii = 0; ii < N_cell; ++ii) {
        bx[ii] = REAL(_B)[ii];
    }
    
    int i;
    
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < N_feature; ++i) {
        int ii = INTEGER(idx)[i]  -1;
        if (ap[ii] == ap[ii+1]) {
            REAL(LXval)[i] = 0;
            // REAL(LYval)[i] = 0;
            REAL(Rval)[i]  = 0;
            REAL(Dval)[i]  = 0;
            REAL(Tval)[i]  = 0;
            continue;
        }
        
        double *tmpa = R_Calloc(N_cell, double);
        double *tmpb = R_Calloc(N_cell, double);
        memset(tmpa, 0, sizeof(double)*N_cell);
        memset(tmpb, 0, sizeof(double)*N_cell);
        double mna = 0,
            mnb = 0;
        int j;
        for (j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            int cid = ai[j];
            tmpa[cid] = ax[j];
        }
        
        for (j = 0; j < N_cell; ++j) tmpb[j] = bx[j];

        for (j = 0; j < N_cell; ++j) {
            if (sens) {
                tmpb[j] = tmpb[j] - tmpa[j];
                if (tmpb[j] < 0) tmpb[j] = 0;
            }
            tmpb[j] = tmpb[j]/REAL(cs)[j];
            mnb += tmpb[j];
        }
        mnb = mnb/N_cell;

        for (j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;            
            int cid = ai[j];
            tmpa[cid] = ax[j]/REAL(cs)[cid]*scale_factor;
            mna += tmpa[cid];
        }
        mna = mna/N_cell;
        
        double *tmpa_s = R_Calloc(N_cell, double);
        double *tmpb_s = R_Calloc(N_cell, double);
        smooth_W(tmpa, tmpa_s, N_cell, W);
        smooth_W(tmpb, tmpb_s, N_cell, W);
        double mna_s = 0,
            mnb_s = 0;
        for (j = 0; j < N_cell; ++j) {
            mna_s += tmpa_s[j];
            mnb_s += tmpb_s[j];
        }
        mna_s = mna_s/(double)N_cell;
        mnb_s = mnb_s/(double)N_cell;
        
        double Lx1 = 0,
            Lx2 = 0,
            //Ly1 = 0,
            //Ly2 = 0,
            ra  = 0,
            rb1 = 0,
            rb2 = 0;
        for (j = 0; j < N_cell; ++j) {
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
        // double Ly = Ly1/Ly2;
        double r = ra/(rb1*rb2);
        double e = sqrt(Lx) * (1-r);
#pragma omp critical
        {
            REAL(LXval)[i] = Lx;
            //REAL(LYval)[i] = Ly;
            REAL(Rval)[i]  = r;
            REAL(Dval)[i]  = e;
        }

        // mean, var
        double mean = 0,
            var = 0;
        
        double *es = R_Calloc(perm, double);
        int k;
        for (k = 0; k < perm; ++k) {
            shuffle(tmpa, ris[k], N_cell);
            smooth_W(tmpa, tmpa_s, N_cell, W);
            mna_s = 0;
            for (j = 0; j < N_cell; ++j) {
                mna_s += tmpa_s[j];
            }
            mna_s = mna_s/(double)N_cell;
            Lx1 = 0;
            ra = 0;
            rb1 = 0;
            for (j = 0; j < N_cell; ++j) {
                Lx1 += pow(tmpa_s[j]-mna,2);
                tmpa_s[j] = tmpa_s[j] - mna_s;
                ra += tmpa_s[j] * tmpb_s[j];
                rb1 += pow(tmpa_s[j],2);
            }
            rb1 = sqrt(rb1);
            Lx = Lx1/Lx2;
            r = ra/(rb1*rb2);
            es[k] = sqrt(Lx) *(1-r);
            mean += es[k];
        }
        mean = mean/perm;

        for (k = 0; k < perm; ++k) {
            var += pow((es[k] -mean),2);
        }
        var = sqrt(var/perm);

        double t = (e - mean)/var;
        
        R_Free(es);
        R_Free(tmpa_s);
        R_Free(tmpb_s);
        R_Free(tmpa);
        R_Free(tmpb);
        
#pragma omp critical
        {
            REAL(Tval)[i]  = t;
            REAL(Mval)[i]  = mean;
            REAL(Vval)[i]  = var;
        }
}

    for (int pi = 0; pi < perm; ++pi) R_Free(ris[pi]);
    R_Free(ris);
    R_Free(bx);
    
    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&W, &c);

    SEXP ta = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ta, 0, LXval);
    //SET_VECTOR_ELT(ta, 1, LYval);
    SET_VECTOR_ELT(ta, 1, Rval);
    SET_VECTOR_ELT(ta, 2, Dval);
    SET_VECTOR_ELT(ta, 3, Tval);
    SET_VECTOR_ELT(ta, 4, Mval);
    SET_VECTOR_ELT(ta, 5, Vval);

    UNPROTECT(7);
    return ta;
}
double *sp_colsums(CHM_SP A, Rboolean mn)
{
    double *cs = R_Calloc(A->ncol, double);
    const int * ap = (int*)A->p;
    //const int * ai = (int*)A->i;
    const double * ax = (double*)A->x;

    for (int i = 0; i < A->ncol; ++i) {
        if (ap[i] == ap[i+1]) {
            cs[i] = 0;
            continue;
        }
        for (int j = ap[i]; j < ap[i+1]; ++j) {
            if (ISNAN(ax[i])) continue;
            cs[i] += ax[j];
        }

        if (mn) cs[i] = cs[i]/A->ncol;
    }
    return cs;
}
/*
 * workspace:
 *  tmpa (A->nrow) scaled expression for cells
 *  allocate temporary copy of W' and W + W'
 *
*/
SEXP autocorrelation_test(SEXP _A, SEXP _W, SEXP _random, SEXP _threads)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP W = AS_CHM_SP__(_W);
    Rboolean rand = asLogical(_random);
    int n_thread = asInteger(_threads);
    
    if (A->stype) return mkString("A cannot be symmetric");
    if (W->stype) return mkString("W cannot be symmetric");
    
    double one[] = {1, 0};

    R_CheckStack();

    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (A->ncol < 2) return mkString("Too few cells."); // to do

    A = M_cholmod_transpose(A, (int)A->xtype, &c);
    W = M_cholmod_transpose(W, (int)W->xtype, &c);
    
    int N = A->nrow; // cells
    int nc = A->ncol; // features

    double N1 = N -1;
    double N2 = N -2;
    double N3 = N -3;
    //double N4 = N -4;
    double NN = pow(N, 2);
    double tmp1 = NN - 3 *N + 3;
    double tmp2 = N1 *N2 * N3;
    double tmp3 = NN - N;

    double EI = -1/N1;
    double EI2 = pow(EI,2);
    
    SEXP Ival = PROTECT(allocVector(REALSXP, nc));
    SEXP Cval = PROTECT(allocVector(REALSXP, nc));
    SEXP IXval = PROTECT(allocVector(REALSXP, nc));
    SEXP CXval = PROTECT(allocVector(REALSXP, nc));    
    
    const int *const ap = (int*)A->p;
    const int *const ai = (int*)A->i;
    const double *const ax = (double*)A->x;
    
    const int *const bp = (int*)W->p;
    const int *const bi = (int*)W->i;
    const double *const bx = (double*)W->x;

    double S0 = 0,
        S1 = 0,
        S2 = 0,
        S02;
    
    int i;
    for (i = 0; i < W->nzmax; ++i) {
        S0 += bx[i];
    }
    S02 = pow(S0, 2);
    
    CHM_SP Wt;
    Wt = M_cholmod_transpose(W, (int)W->xtype, &c);

    double *cs = sp_colsums(W, FALSE);
    double *rs = sp_colsums(Wt, FALSE);
    
    CHM_SP WWt;
    WWt = M_cholmod_add(W, Wt, one, one, TRUE, TRUE, &c);
    M_cholmod_free_sparse(&Wt, &c);
    double *WWtx = (double*)WWt->x;
    for (i = 0; i < WWt->nzmax; ++i) S1 += pow(WWtx[i],2);
    S1 = S1/2;

    for (i = 0; i < W->ncol; ++i) S2 += pow((cs[i]+rs[i]),2);
    
    R_Free(cs);
    R_Free(rs);    
    M_cholmod_free_sparse(&WWt, &c);

    const double varI_norm = (NN*S1-N*S2+3*S02)/(S02*(NN-1)) - EI2;
    const double varC_norm = ((2*S1+S2)*N1 - 4 * S02)/(2*(N+1)*S02);

    int ci;
#pragma omp parallel for num_threads(n_thread)
    for (ci = 0; ci < nc; ++ci) { // for each feature
        double *tmpa = R_Calloc(N, double);
        // progress(nc);
        
        if (ap[ci] == ap[ci+1]) {
            REAL(Ival)[ci] = 0;
            continue;   
        }
        
        // memset(tmpa, 0, sizeof(double)*N);

        int j;
        double varI = 0,
            varC = 0,
            mn = 0,
            xix = 0,
            xij = 0,
            m2 = 0,
            m4 = 0,
            b2,
            sigma2 = 0,
            C = 0;
        
        for (j = ap[ci]; j < ap[ci+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            mn += ax[j];
        }
        mn = mn/(double)N;
        
        for (j = ap[ci]; j < ap[ci+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            tmpa[ai[j]] = ax[j];
        }
        
        // for (j = 0; j < N; ++j) tmpa[j] = tmpa[j] -mn;        
        for (j = 0; j < N; ++j) { // for cells
            xix += pow(tmpa[j]-mn,2);
            m4 += pow(tmpa[j]-mn, 4);
        }

        sigma2 = xix/N1;
        
        int wi;
        for (wi = 0; wi < N; ++wi) {
            if (bp[wi] == bp[wi+1]) continue;
            int wj;
            for (wj = bp[wi]; wj < bp[wi+1]; ++wj) {
                // Rprintf("W : %f, Zi : %f, Zj : %f\n", bx[wj], tmpa[wi], tmpa[bi[wj]]);
                C += bx[wj] * pow(tmpa[wi] - tmpa[bi[wj]],2);
            }
        }

        if (rand) { // random
            m2 = xix/N;
            m4 = m4/N;
            b2 = m4/pow(m2,2);
            varI = N*(tmp1*S1 - N*S2 + 3*S02) - b2*(tmp3*S1 - 2*N*S2 + 6*S02);
            varI = varI/(tmp2*S02);
            varI = varI - EI2;
            varC = N1*S1*(NN-3*N+3-N1*b2) - N1*S2*(NN+3*N-6-(NN-N+2)*b2)/4 + S02*(NN-3-N1*N1*b2);
            varC = varC/(N*N2*N3*S02);
        }
        
        for (int i = 0; i < N; ++i) {
            for (j = bp[i]; j < bp[i+1]; ++j) {
                if (ISNAN(bx[j])) continue;
                int cid = bi[j];
                if (i == cid) continue; // i != j
                xij += (bx[j] * (tmpa[i]-mn) * (tmpa[cid]-mn));
            }
        }

        double I = N*xij/(xix*S0);
        //Rprintf("wC : %f, sigma : %f, S0 : %f\n", C, sigma2, S0);
        C = C/(2*sigma2*S0);
        //Rprintf("Cvar : %f, %f\n", varC, varC_norm);
        R_Free(tmpa);
#pragma omp critical
        {
            REAL(Ival)[ci] = I;
            REAL(Cval)[ci] = C;
            REAL(IXval)[ci] = rand ? (I - EI)/sqrt(varI) : (I - EI)/sqrt(varI_norm);
            REAL(CXval)[ci] = rand ? (C - 1)/sqrt(varC) : (C - 1)/sqrt(varC_norm);
        }
    }


    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&W, &c);
    
    SEXP ta = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ta, 0, Ival);
    SET_VECTOR_ELT(ta, 1, Cval);
    SET_VECTOR_ELT(ta, 2, IXval);
    SET_VECTOR_ELT(ta, 3, CXval);

    UNPROTECT(5);
    return ta;
}

