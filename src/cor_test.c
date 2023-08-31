#include<R.h>
#include<Rdefines.h>
#include<Rinternals.h>
#include "Matrix_stubs.c"

#include <omp.h>

void progress(int p) {
    static int displayed = -1;
    static int count = 0;
    static char bar[] = "================================================== ";

    count++;

    int r = count*100/p;
    r/=2;
    
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
cholmod_common c;

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
        double sd = 0.0;

        for (j = ap[ci]; j < ap[ci+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            mn += ax[j];
            tmpa[ai[j]] = ax[j];
        }
        mn = mn/N;

        for (j = 0; j < N; ++j) tmpa[j] = tmpa[j] - mn;
        
        mn = 0;
        sd = 0;        
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

void smooth_W(double * const a, double *s, const int N, const_CHM_SP W)
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
    
    double one[] = {1, 0};

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

    const int const *ap = (int*)A->p;
    const int const *ai = (int*)A->i;
    const double const *ax = (double*)A->x;

    const int const *wp = (int*)W->p;
    const int const *wi = (int*)W->i;
    const double const *wx = (double*)W->x;

    int *hits = R_Calloc(N_feature, int);
    
    int *steps = R_Calloc(N_feature, int);
    int N = N_feature;

    int n_step = perm/CACHE_PER_BATCH+1;
    
    for (int i = 0; i < N; ++i) steps[i] = i;
    
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
                sd = 0,
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
            
            int k;
            for (k = 0; k < N_cell; ++k) {
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
                
                int pi;
                for (pi = 0; pi < CACHE_PER_BATCH; ++pi) {
                    shuffle(tmp, ris[pi], N_cell);
                    int k;
                    for (k = 0; k < N_cell; ++k) {
                        for (j = wp[k]; j < wp[k+1]; ++j) {
                            if (ISNAN(wx[j])) continue;
                            int cid = wi[j];
                            if (k == cid) continue; // i != j
                            xij += wx[j] * tmp[k] * tmp[cid];
                        }
                    }
                    Is[pi] = xij/xix;
                }
                for (k = 0; k < CACHE_PER_BATCH; ++k) {
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
            int pi;
            for (pi = 0; pi < CACHE_PER_BATCH; ++pi) R_Free(ris[pi]);
            R_Free(ris);
        }

        int j = 0;
        for (i = 0; i < N_feature; ++i) {
            if (hits[i] > MIN_HIT) {
                if (REAL(Pval)[i] == 0) {
                    REAL(Pval)[i] == (double)hits[i]/((step+1)*CACHE_PER_BATCH);
                }
            } else {
                steps[j++] = i;
            }
        }
        N = j;
    }

    for (int i = 0; i < N; ++i) {
        int idx = steps[i];
        REAL(Pval)[idx] == (double)hits[idx]/(n_step*CACHE_PER_BATCH);
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
 
SEXP E_test(SEXP _A, SEXP _B, SEXP _W, SEXP _permut, SEXP _threads)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP B = AS_CHM_SP__(_B);
    CHM_SP W = AS_CHM_SP__(_W);
    int perm = asInteger(_permut);
    int n_thread = asInteger(_threads);
    
    if (A->stype) return mkString("A cannot be symmetric");
    if (B->stype) return mkString("B cannot be symmetric");
    if (W->stype) return mkString("W cannot be symmetric");
    
    double one[] = {1, 0};

    R_CheckStack();

    if (A->ncol != B->ncol || A->nrow != B->nrow) return mkString("A and B do not match");
    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (A->ncol < 2) return mkString("Too few cells."); // to do

    A = M_cholmod_transpose(A, (int)A->xtype, &c);
    B = M_cholmod_transpose(B, (int)B->xtype, &c);
    W = M_cholmod_transpose(W, (int)W->xtype, &c);
    
    const int N_cell = A->nrow;
    const int N_feature = A->ncol;

    SEXP LXval = PROTECT(allocVector(REALSXP, N_feature));
    SEXP LYval = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Rval = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Eval = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Pval = PROTECT(allocVector(REALSXP, N_feature));
    
    const int * const ap = (int*)A->p;
    const int * const ai = (int*)A->i;
    const double * const ax = (double*)A->x;

    const int *const bp = (int*)B->p;
    const int *const bi = (int*)B->i;
    const double *const bx = (double*)B->x;
    
    int *hits = R_Calloc(N_feature, int);
    for (int i = 0; i < N_feature; ++i) hits[i] = 0;
    
    int n_step = perm/CACHE_PER_BATCH;
    if (n_step < 1) n_step = 1;
    
    int *steps = R_Calloc(N_feature, int);
    for (int i = 0; i < N_feature; ++i) steps[i] = i;
    
    int N = N_feature;
    
    for (int step = 0; step < n_step; ++step) {
        int **ris = NULL;
        if (perm > 1) {
            REprintf("Generate weight matrix idx ... ");
            ris = R_Calloc(CACHE_PER_BATCH,int*);
            int pi;
            for (pi = 0; pi < CACHE_PER_BATCH; ++pi) {
                ris[pi] = random_idx(N_cell);
            }
            REprintf("Finished.\n");
        }
    
        int i;
#pragma omp parallel for num_threads(n_thread)
        for (i = 0; i < N; ++i) {
            int idx = steps[i];
            if (ap[idx] == ap[idx+1] || bp[idx] == bp[idx+1]) {
                REAL(LXval)[i] = 0;
                REAL(LYval)[i] = 0;
                REAL(Rval)[i] = 0;
                REAL(Eval)[i] = 0;
                REAL(Pval)[i] = 0;
                continue;
            }

            double *tmpa = R_Calloc(N_cell, double);
            double *tmpb = R_Calloc(N_cell, double);

            memset(tmpa, 0, sizeof(double)*N_cell);
            memset(tmpb, 0, sizeof(double)*N_cell);
            
            int j;
            double mna = 0,
                mnb = 0;
            
            for (j = ap[idx]; j < ap[idx+1]; ++j) {
                if (ISNAN(ax[j])) continue;
                tmpa[ai[j]] = ax[j];
                mna += ax[j];
            }
            mna = mna/N_cell;
            
            for (j = bp[idx]; j < bp[idx+1]; ++j) {
                if (ISNAN(bx[j])) continue;
                tmpb[bi[j]] = bx[j];
                mnb += bx[j];
            }
            mnb = mnb/N_cell;

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
                ra = 0,
                rb1 = 0,
                rb2 = 0;
            
            for (j = 0; j < N_cell; ++j) {
                Lx1 += pow(tmpa_s[j]-mna,2);
                Lx2 += pow(tmpa[j]-mna,2);
                Ly1 += pow(tmpb_s[j]-mnb,2);
                Ly2 += pow(tmpb[j]-mnb,2);
                
                tmpa_s[j] = tmpa_s[j] - mna_s;
                tmpb_s[j] = tmpb_s[j] - mnb_s;
                
                ra += tmpa_s[j] * tmpb_s[j];
                rb1 += pow(tmpa_s[j],2);
                rb2 += pow(tmpb_s[j],2);
            }
            
            rb1 = sqrt(rb1);
            rb2 = sqrt(rb2);
            
            double Lx = Lx1/Lx2;
            double Ly = Ly1/Ly2;
            double r = ra/(rb1*rb2);
            double e = sqrt(Lx) * (1-r);
            
            if (perm > 1) {
                double *es = R_Calloc(CACHE_PER_BATCH, double);
                int k;
                for (k = 0; k < CACHE_PER_BATCH; ++k) {
                    //shuffle_double_arr(tmpa, N_cell);
                    shuffle(tmpa, ris[k], N_cell);
                    //shuffle(tmpb, ris[k], N_cell);
                    
                    smooth_W(tmpa, tmpa_s, N_cell, W);
                    //smooth_W(tmpb, tmpb_s, N_cell, W);
                    
                    mna_s = 0;
                    //mnb_s = 0;
                    
                    for (j = 0; j < N_cell; ++j) {
                        mna_s += tmpa_s[j];
                        //mnb_s += tmpb_s[j];
                    }
                    mna_s = mna_s/(double)N_cell;
                    //mnb_s = mnb_s/(double)N_cell;
                    
                    Lx1 = 0;
                    //Ly1 = 0;
                    ra = 0;
                    rb1 = 0;
                    //rb2 = 0;
                    
                    for (j = 0; j < N_cell; ++j) {
                        Lx1 += pow(tmpa_s[j]-mna,2);
                        // Ly1 += pow(tmpb_s[j]-mnb,2);
                        tmpa_s[j] = tmpa_s[j] - mna_s;
                        // tmpb_s[j] = tmpb_s[j] - mnb_s;
                        
                        ra += tmpa_s[j] * tmpb_s[j];
                        rb1 += pow(tmpa_s[j],2);
                        //rb2 += pow(tmpb_s[j],2);
                    }
                    
                    rb1 = sqrt(rb1);
                    // rb2 = sqrt(rb2);
                    
                    Lx = Lx1/Lx2;
                    //Ly = Ly1/Ly2;
                    
                    r = ra/(rb1*rb2);
                    //es[k] = sqrt(Lx) * sqrt(Ly) *(1-r);
                    es[k] = sqrt(Lx) *(1-r);
                    /* Rprintf("Lx1 : %f %f, Lx2 : %f %f, ra : %f %f, r : %f %f, mn : %f %f\n", */
                    /*         Lx1, Lx1_p, Lx2, Lx2_p, ra, ra_p, r , r_p, mna, mna_p); */
                }
                
                for (k = 0; k < CACHE_PER_BATCH; ++k) {
                    //Rprintf("%f\t%f\n",e, es[k]);
                    if (es[k]>=e) hits[idx] +=1;
                }
                R_Free(es);
            }

            R_Free(tmpa_s);
            R_Free(tmpb_s);
            R_Free(tmpa);
            R_Free(tmpb);
            
#pragma omp critical
            if (step == 0) {
                REAL(LXval)[i] = Lx;
                REAL(LYval)[i] = Ly;
                REAL(Rval)[i] = r;
                REAL(Eval)[i] = e;
                // REAL(Pval)[i] = p;
            }
        }
        
        if (perm > 1) {
            int pi;
            for (pi = 0; pi < CACHE_PER_BATCH; ++pi) R_Free(ris[pi]);
            R_Free(ris);
        }

        int j = 0;
        for (i = 0; i < N_feature; ++i) {
            if (hits[i] == -1) continue;
            
            if (hits[i] > MIN_HIT) {
                REAL(Pval)[i] == (double)hits[i];///((step+1)*CACHE_PER_BATCH);
                hits[i] = -1;
            } else {
                steps[j++] = i;
            }
        }
        N = j;

        Rprintf("N : %d\n", N);
        if (N == 0) break;
    }

    for (int i = 0; i < N; ++i) {
        int idx = steps[i];
        REAL(Pval)[idx] == (double)hits[idx];///(n_step*CACHE_PER_BATCH);
    }
    R_Free(hits);
    R_Free(steps);
    
    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&B, &c);
    M_cholmod_free_sparse(&W, &c);

    SEXP ta = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ta, 0, LXval);
    SET_VECTOR_ELT(ta, 1, LYval);
    SET_VECTOR_ELT(ta, 2, Rval);
    SET_VECTOR_ELT(ta, 3, Eval);
    SET_VECTOR_ELT(ta, 4, Pval);

    UNPROTECT(6);
    return ta;
}
double *sp_colsums(CHM_SP A, Rboolean mn)
{
    double *cs = R_Calloc(A->ncol, double);
    int *ap = (int*)A->p;
    int *ai = (int*)A->i;
    double *ax = (double*)A->x;
    int i;
    for (i = 0; i < A->ncol; ++i) {
        if (ap[i] == ap[i+1]) {
            cs[i] = 0;
            continue;
        }
        int j;
        for (j = ap[i]; j < ap[i+1]; ++j) {
            if (ISNAN(ax[i])) continue;
            cs[i] += ax[j];
        }

        if (mn) cs[i]/A->ncol;
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
    double N4 = N -4;
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
    
    int *ap = (int*)A->p;
    int *ai = (int*)A->i;
    double *ax = (double*)A->x;

    int *bp = (int*)W->p;
    int *bi = (int*)W->i;
    double *bx = (double*)W->x;

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

    double varI_norm = (NN*S1-N*S2+3*S02)/(S02*(NN-1)) - EI2;
    double varC_norm = ((2*S1+S2)*N1 - 4 * S02)/(2*(N+1)*S02);

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
        
        int i;
        for (i = 0; i < N; ++i) {
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

