#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
//#include<Rinternals.h>

#ifndef Matrix_stubs
#define Matrix_stubs
#include "Matrix_stubs.c"
#endif

#include <assert.h>

// # nocov end
static cholmod_common c;

#define CACHE_PER_BATCH 10000
#define MIN_HIT 1

SEXP openmp_support()
{
#ifdef _OPENMP
    return ScalarLogical(1);
#endif
    return ScalarLogical(0);
}
double * shuffle_double_arr(double *a, const int n)
{
    double *s = R_Calloc(n, double);
    int i;
    for (i = 0; i < n-1; ++i) {
        int j = i + rand() / (RAND_MAX / (n - i) + 1);
        s[i] = a[j];
        s[j] = a[i];
    }
    return(s);
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
void shuffle(double tmp[], int const idx[], const int n)
{
    double aux[n];
    int i;
    for (i = 0; i < n; i++) aux[idx[i]] = tmp[i];
    for (i = 0; i < n; i++) tmp[i] = aux[i];
}
/*
SEXP cor_test(SEXP _A, SEXP _B, SEXP _trans)
{
    Rboolean tr = asLogical(_trans);
    
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP B = AS_CHM_SP__(_B);

    CHM_SP A2 = NULL;
    CHM_SP B2 = NULL;

    int freeA = 0;
    int freeB = 0;
    
    if (A->ncol != B->ncol || A->nrow != B->nrow) return mkString("Unequal matrix size.");
    if (A->ncol < 2) return mkString("Too few cells.");

    if (A->stype) {
        A2 = M_cholmod_copy(A, 0, TRUE, &c);
        if (c.status < CHOLMOD_OK) return mkString("Out of memory!");
        A = A2;
        freeA = 1;
    }

    if (B->stype) {
        B2 = M_cholmod_copy(B, 0, TRUE, &c);
        if (c.status < CHOLMOD_OK) return mkString("Out of memory!");        
        B = B2;
        freeB = 1;
    }
    if (tr) {
        A2 = M_cholmod_transpose(A, (int)A->xtype, &c);
        B2 = M_cholmod_transpose(B, (int)B->xtype, &c);
        if (freeA) M_cholmod_free_sparse(&A, &c);
        if (freeB) M_cholmod_free_sparse(&B, &c);

        A = A2;
        B = B2;
        freeA = 1;
        freeB = 1;
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

    if (freeA) M_cholmod_free_sparse(&A, &c);
    if (freeB) M_cholmod_free_sparse(&B, &c);

    UNPROTECT(1);
    return rval;
}
*/
void smooth_W(double * const a, double *s, const int N, CHM_SP W)
{
    int *wp = (int*)W->p;
    int *wi = (int*)W->i;
    double *wx = (double*)W->x;

    memset(s, 0, N*sizeof(double));

    int i;
    for (i = 0; i < N; ++i) {
        if (wp[i] == wp[i+1]) continue;
        //double test = 0;
        int j;
        for (j = wp[i]; j < wp[i+1]; ++j) {
            int idx = wi[j];
            s[i] += wx[j] * a[idx];
            //test += wx[j];
        }
        //Rprintf("W: %f", test);
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
    
    int n_cell = A->nrow;
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
                ris[pi] = random_idx(n_cell);
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
            
            double *tmp = R_Calloc(n_cell, double);
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
            mn = mn/n_cell;
            
            for (j = 0; j < n_cell; ++j) {
                tmp[j] = tmp[j] - mn;
                xix += pow(tmp[j], 2);
            }
            
            for (int k = 0; k < n_cell; ++k) {
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
                    shuffle(tmp, ris[pi], n_cell);
                    for (int k = 0; k < n_cell; ++k) {
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
/*
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
    const int n_cell = A->nrow;
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
        
        double *tmpa = R_Calloc(n_cell, double);
        memset(tmpa, 0, sizeof(double)*n_cell);
        int j;
        for (j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            int cid = ai[j];
            tmpa[cid] = log(ax[j]/REAL(cs)[cid]*scale_factor + 1);
        }
        
        double *tmpa_s = R_Calloc(n_cell, double);
        smooth_W(tmpa, tmpa_s, n_cell, W);
        SEXP val = PROTECT(allocVector(REALSXP, n_cell));
        for (j = 0; j < n_cell; ++j) {
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

SEXP D_score_lite(SEXP _A, SEXP _B, SEXP _W)
{
    int n_cell = length(_A);
    assert(n_cell == length(_B));
    CHM_SP W = AS_CHM_SP__(_W);
    if (W->stype) return mkString("W cannot be symmetric");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");    
    if (W->ncol != n_cell) return mkString("Inconsistant cell number.");

    double *tmpa = R_Calloc(n_cell, double);
    double *tmpb = R_Calloc(n_cell, double);
    memset(tmpa, 0, sizeof(double)*n_cell);
    memset(tmpb, 0, sizeof(double)*n_cell);
    double mna = 0, mnb = 0;

    for (int i = 0; i < n_cell; ++i) {
        tmpa[i] = REAL(_A)[i];
        tmpb[i] = REAL(_B)[i];
        mna = mna + tmpa[i];
        mnb = mnb + tmpb[i];
    }
    mna = mna/n_cell;
    mnb = mnb/n_cell;

    double *tmpa_s = R_Calloc(n_cell, double);
    double *tmpb_s = R_Calloc(n_cell, double);
    smooth_W(tmpa, tmpa_s, n_cell, W);
    smooth_W(tmpb, tmpb_s, n_cell, W);
    double mna_s = 0, mnb_s = 0;
        
    for (int j = 0; j < n_cell; ++j) {
        mna_s += tmpa_s[j];
        mnb_s += tmpb_s[j];
    }
    mna_s = mna_s/(double)n_cell;
    mnb_s = mnb_s/(double)n_cell;

    double Lx1 = 0,
        Lx2 = 0,
        Ly1 = 0,
        Ly2 = 0,
        ra  = 0,
        rb1 = 0,
        rb2 = 0;
    for (int j = 0; j < n_cell; ++j) {
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
    
    R_Free(tmpa);
    R_Free(tmpb);
    R_Free(tmpa_s);
    R_Free(tmpb_s);

    rb1 = sqrt(rb1);
    rb2 = sqrt(rb2);
    
    double Lx = Lx1/Lx2;
    double r = ra/(rb1*rb2);
    
    SEXP result = PROTECT(allocVector(REALSXP, 1));
    REAL(result)[0] = sqrt(Lx) * (1-r);
    UNPROTECT(1);
    
    return result; 
}
*/
SEXP D_test(SEXP _A,
            SEXP _B,
            SEXP _W,
            SEXP _method,
            SEXP _permut,
            SEXP _threads,
            SEXP idx,
            SEXP bidx,
            SEXP cs,
            SEXP _factor,
            SEXP _mode,
            SEXP _scale,
            SEXP _norm,
            SEXP _seed,
            SEXP _debug)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP B = AS_CHM_SP__(_B);
    CHM_SP W = AS_CHM_SP__(_W);

    CHM_SP A2;
    CHM_SP B2;
    
    int freeA = 0;
    int freeB = 0;
    
    const int perm = asInteger(_permut);
    int n_thread = asInteger(_threads);
    const int method = asInteger(_method);
    
    const int scale_factor = asInteger(_factor);
    Rboolean scale = asLogical(_scale);
    Rboolean norm = asLogical(_norm);
    Rboolean debug = asLogical(_debug);
    
    const int seed = asInteger(_seed);
    srand(seed);
    
    int mode = asInteger(_mode);
    if (mode != 1 && mode != 2 && mode != 3) return mkString("Unsupported mode!");

    if (debug) {
        Rprintf("mode %d\n", mode);
    }
    if (W->stype) return mkString("W cannot be symmetric");
    if (A->stype) {
        A2 = M_cholmod_copy(A, 0, TRUE, &c);
        if (c.status < CHOLMOD_OK) return mkString("Out of memory!");
        A = A2;
        freeA = 1;
    }
    if (B->stype) {
        B2 = M_cholmod_copy(B, 0, TRUE, &c);
        if (c.status < CHOLMOD_OK) return mkString("Out of memory!");        
        B = B2;
        freeB = 1;
    }

    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (A->ncol < 2) return mkString("Too few cells."); // to do

    A2 = M_cholmod_transpose(A, (int)A->xtype, &c);
    B2 = M_cholmod_transpose(B, (int)B->xtype, &c);

    if (freeA) M_cholmod_free_sparse(&A, &c);
    if (freeB) M_cholmod_free_sparse(&B, &c);
    A = A2;
    B = B2;
    
    const int n_cell = A->nrow;
    const int N_feature = length(idx);

    if (debug) {
        Rprintf("n_cell, %d, n_feature, %d\n", n_cell, N_feature);
        n_thread = 1; // disable multithreads
    }
    assert (length(bidx) == N_feature);
    
    // SEXP LXval = PROTECT(allocVector(REALSXP, N_feature));
    // SEXP LYval = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Rval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Dval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Tval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Mval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Vval  = PROTECT(allocVector(REALSXP, N_feature));
    
    const int * ap = (int*)A->p;
    const int * ai = (int*)A->i;
    const double * ax = (double*)A->x;

    const int * bp = (int*)B->p;
    const int * bi = (int*)B->i;
    const double * bx = (double*)B->x;
    
    int **ris = NULL;
    ris = R_Calloc(perm,int*);
    for (int pi = 0; pi < perm; ++pi) {
        ris[pi] = random_idx(n_cell);
    }
    
    R_CheckUserInterrupt();

    int i;
    
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < N_feature; ++i) {
        int ii = INTEGER(idx)[i]  -1;
        int ij = INTEGER(bidx)[i] -1;
        if (ap[ii] == ap[ii+1] || bp[ij] == bp[ij+1]) {
            REAL(Rval)[i]  = NA_REAL;
            REAL(Dval)[i]  = NA_REAL;
            REAL(Tval)[i]  = NA_REAL;
            REAL(Mval)[i]  = NA_REAL;
            REAL(Vval)[i]  = NA_REAL;
            continue;
        }

        double *tmpa = R_Calloc(n_cell,double);
        double *tmpb = R_Calloc(n_cell,double);
        
        // double *tmpa = alloca(n_cell *sizeof(double));
        // double *tmpb = alloca(n_cell *sizeof(double));
        // R_CheckStack();
        
        memset(tmpa, 0, sizeof(double)*n_cell);
        memset(tmpb, 0, sizeof(double)*n_cell);
        
        double mna = 0, mnb = 0;
        for (int j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            int cid = ai[j];
            tmpa[cid] = ax[j];
        }

        for (int j = bp[ij]; j < bp[ij+1]; ++j) {
            if (ISNAN(bx[j])) continue;
            int cid = bi[j];
            tmpb[cid] = bx[j];
        }

        for (int j = 0; j < n_cell; ++j) {
            if (mode == 2) {
                tmpb[j] = tmpb[j] - tmpa[j];
                if (tmpb[j] < 0) tmpb[j] = 0;
            } else if (mode == 3) {
                tmpb[j] = tmpb[j] + tmpa[j];
            }

            if (norm) {
                tmpb[j] = log(tmpb[j]/REAL(cs)[j]*scale_factor + 1);
            }
            mnb += tmpb[j];
        }
        mnb = mnb/n_cell;

        for (int j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;            
            int cid = ai[j];            
            if (norm) {
                tmpa[cid] = log(ax[j]/REAL(cs)[cid]*scale_factor + 1);
            }
            mna += tmpa[cid];
        }
        mna = mna/n_cell;

        if (debug) {
            Rprintf("mean_a, %f, mean_b, %f\n", mna, mnb);
        }
        if (scale) {
            double sd1 = 0;
            double sd2 = 0;
            for (int j = 0; j < n_cell; ++j) {
                sd1 += pow(tmpa[j] - mna, 2);
                sd2 += pow(tmpb[j] - mnb, 2);
            }
            sd1 = sqrt(sd1/n_cell);
            sd2 = sqrt(sd2/n_cell);

            for (int j = 0; j < n_cell; ++j) {
                tmpa[j] = (tmpa[j] - mna)/sd1;
                tmpb[j] = (tmpb[j] - mnb)/sd2;
            }

            mna = 0;
            mnb = 0;
        }
        
        double *tmpa_s = R_Calloc(n_cell, double);
        double *tmpb_s = R_Calloc(n_cell, double);
        // double *tmpa_s = alloca(n_cell *sizeof(double));
        // double *tmpb_s = alloca(n_cell *sizeof(double));
        smooth_W(tmpa, tmpa_s, n_cell, W);
        smooth_W(tmpb, tmpb_s, n_cell, W);
        double mna_s = 0;
        double mnb_s = 0;
        for (int j = 0; j < n_cell; ++j) {
            mna_s += tmpa_s[j];
            mnb_s += tmpb_s[j];
        }
        mna_s = mna_s/n_cell;
        mnb_s = mnb_s/n_cell;

        if (debug) {
            Rprintf("smooth mean_a, %f, smooth mean_b, %f\n", mna_s, mnb_s);
        }

        double Lx1 = 0, Lx2 = 0, Ly1 = 0, Ly2 = 0, ra = 0, rb1 = 0, rb2 = 0;
        
        for (int j = 0; j < n_cell; ++j) {
            //Lx1 += pow(tmpa_s[j]-mna_s,2);
            Lx1 += pow(tmpa_s[j]-mna,2);
            Lx2 += pow(tmpa[j]-mna,2);

            if (method != 1) {
                Ly1 += pow(tmpb_s[j]-mnb,2);
                Ly2 += pow(tmpb[j]-mnb,2);
            }

            tmpa_s[j] = tmpa_s[j] - mna_s;
            tmpb_s[j] = tmpb_s[j] - mnb_s;
            ra  += tmpa_s[j] * tmpb_s[j];
            rb1 += pow(tmpa_s[j],2);
            rb2 += pow(tmpb_s[j],2);
        }
        
        rb1 = sqrt(rb1);
        rb2 = sqrt(rb2);
        
        double Lx = Lx1/Lx2;        
        double r = ra/(rb1*rb2);
        double Ly, e;
        if (method == 1) {
            e = sqrt(Lx) * (1-r);
        } else if (method == 2) {
            Ly = Ly1/Ly2;
            e = sqrt(Lx) * sqrt(Ly) * (1-r);
        } else {
            Ly = Ly1/Ly2;
            e = sqrt(Lx) * sqrt(Ly) * r;
        }
        
        if (debug) {
            Rprintf("Lx, %f, D, %f, r, %f\n", Lx, e, r);
        }

#pragma omp critical
        {
            // REAL(LXval)[i] = Lx;
            // REAL(LYval)[i] = Ly;
            REAL(Rval)[i]  = r;
            REAL(Dval)[i]  = e;
        }

        // mean, var
        double mean = 0, var = 0;
        double *es = R_Calloc(perm, double);
        // double *es = alloca(perm *sizeof(double));
        // R_CheckStack();
        
        for (int k = 0; k < perm; ++k) {
            shuffle(tmpa, ris[k], n_cell);
            smooth_W(tmpa, tmpa_s, n_cell, W);
            mna_s = 0;
            for (int j = 0; j < n_cell; ++j) {
                mna_s += tmpa_s[j];
            }
            mna_s = mna_s/n_cell;
            Lx1 = 0;
            ra = 0;
            rb1 = 0;
            for (int j = 0; j < n_cell; ++j) {
                // Lx1 += pow(tmpa_s[j]-mna_s,2);
                Lx1 += pow(tmpa_s[j]-mna, 2);
                tmpa_s[j] = tmpa_s[j] - mna_s;
                rb1 += pow(tmpa_s[j], 2);
            }
            rb1 = sqrt(rb1);
            Lx = Lx1/Lx2;

            if (method != 1) {
                shuffle(tmpb, ris[k], n_cell);
                smooth_W(tmpb, tmpb_s, n_cell, W);
                mnb_s = 0;
                for (int j = 0; j < n_cell; ++j) {
                    mnb_s += tmpb_s[j];
                }
                mnb_s = mnb_s/n_cell;
                Ly1 = 0;
                rb2 = 0;
                for (int j = 0; j < n_cell; ++j) {
                    Ly1 += pow(tmpb_s[j]-mnb, 2);
                    tmpb_s[j] = tmpb_s[j] - mnb_s;
                    rb2 += pow(tmpb_s[j], 2);
                }
                rb2 = sqrt(rb2);
                Ly = Ly1/Ly2;
            }
            
            for (int j = 0; j < n_cell; ++j) {
                ra += tmpa_s[j] * tmpb_s[j];
            }
            
            r = ra/(rb1*rb2);
            if (method == 1) {
                es[k] = sqrt(Lx) *(1-r);
            } else if (method == 2) {
                es[k] = sqrt(Lx) * sqrt(Ly) *(1-r);
            } else {
                es[k] = sqrt(Lx) * sqrt(Ly) * r;
            }

            mean += es[k];
        }
        mean = mean/perm;

        for (int k = 0; k < perm; ++k) {
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

        if (debug) {
            Rprintf("t, %f, mean, %f, var, %f \n", t, mean, var);
        }

    }

    for (int pi = 0; pi < perm; ++pi) R_Free(ris[pi]);
    R_Free(ris);

    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&B, &c);
    //M_cholmod_free_sparse(&W, &c);

    SEXP ta = PROTECT(allocVector(VECSXP, 5));
    // SET_VECTOR_ELT(ta, 0, LXval);
    // SET_VECTOR_ELT(ta, 1, LYval);
    SET_VECTOR_ELT(ta, 0, Rval);
    SET_VECTOR_ELT(ta, 1, Dval);
    SET_VECTOR_ELT(ta, 2, Tval);
    SET_VECTOR_ELT(ta, 3, Mval);
    SET_VECTOR_ELT(ta, 4, Vval);

    UNPROTECT(6);
    return ta;
}
SEXP D_distribution_test(SEXP _A, SEXP _B, SEXP _W, SEXP _permut, SEXP _threads)
{
    int la = length(_A);
    int lb = length(_B);
    
    CHM_SP W = AS_CHM_SP__(_W);
    
    const int perm = asInteger(_permut);
    const int n_thread = asInteger(_threads);
    
    if (W->stype) return mkString("W cannot be symmetric");
    if (la != lb) return mkString("Unequal length of A and B");
   
    const int N = W->nrow;
    
    if (la != N) return mkString("Unequal length of A and W");
    
    SEXP X = PROTECT(allocVector(REALSXP, perm));

    double *a = R_Calloc(N, double);
    double *b = R_Calloc(N, double);
    
    double mn1 = 0;
    double mn2 = 0;
    int i;
    for (i = 0; i < N; ++i) {
        a[i] = REAL(_A)[i];
        b[i] = REAL(_B)[i];
        mn1 += a[i];
        mn2 += b[i];
    }

    mn1 = mn1/N;
    mn2 = mn2/N;

    double *tmpb_s = R_Calloc(N, double);
    smooth_W(b, tmpb_s, N, W);

    double mnb_s = 0;
    for (i = 0; i < N; ++i) {
        mnb_s += tmpb_s[i];
    }
    mnb_s = mnb_s/N;

    double *tmpa_s = R_Calloc(N, double);
    smooth_W(a, tmpa_s, N, W);
    double mna_s = 0;
    int j;
    for (j = 0; j < N; ++j) {
        mna_s += tmpa_s[j];
    }
    mna_s = mna_s/N;
        
    double Lx1 = 0,
        Lx2 = 0,
        ra  = 0,
        rb1 = 0,
        rb2 = 0;
    for (j = 0; j < N; ++j) {
        Lx1 += pow(tmpa_s[j]-mn1,2);
        Lx2 += pow(a[j]-mn1,2);
            
        tmpa_s[j] = tmpa_s[j] - mna_s;
        tmpb_s[j] = tmpb_s[j] - mnb_s;
        ra  += tmpa_s[j] * tmpb_s[j];
        rb1 += pow(tmpa_s[j],2);
        rb2 += pow(tmpb_s[j],2);
    }
        
    rb1 = sqrt(rb1);
    rb2 = sqrt(rb2);

    R_Free(tmpa_s);
    
    double Lx = Lx1/Lx2;
    double r = ra/(rb1*rb2);
    double e = sqrt(Lx) * (1-r);
    Rprintf("r : %f\tD : %f\n", r, e);
        
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < perm; ++i) {
        double *s = shuffle_double_arr(a, N);
        double *tmpa_s = R_Calloc(N, double);
        smooth_W(s, tmpa_s, N, W);
        double mna_s = 0;
        int j;
        for (j = 0; j < N; ++j) {
            mna_s += tmpa_s[j];
        }
        mna_s = mna_s/N;
        
        double Lx1 = 0,
            Lx2 = 0,
            ra  = 0,
            rb1 = 0,
            rb2 = 0;
        for (j = 0; j < N; ++j) {
            Lx1 += pow(tmpa_s[j]-mn1,2);
            Lx2 += pow(s[j]-mn1,2);
            
            tmpa_s[j] = tmpa_s[j] - mna_s;
            tmpb_s[j] = tmpb_s[j] - mnb_s;
            ra  += tmpa_s[j] * tmpb_s[j];
            rb1 += pow(tmpa_s[j],2);
            rb2 += pow(tmpb_s[j],2);
        }
        
        rb1 = sqrt(rb1);
        rb2 = sqrt(rb2);

        R_Free(tmpa_s);
        R_Free(s);
        
        double Lx = Lx1/Lx2;
        double r = ra/(rb1*rb2);
        double e = sqrt(Lx) * (1-r);
#pragma omp critical
        {
            REAL(X)[i]  = e;
        }
    }

    R_Free(tmpb_s);
    R_Free(a);
    R_Free(b);
    UNPROTECT(1);
    return X;
}

/* 
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
    const int n_cell = A->nrow;
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
        ris[pi] = random_idx(n_cell);
    }
    
    double *bx = malloc(n_cell*sizeof(double));
    for (int ii = 0; ii < n_cell; ++ii) {
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
        
        double *tmpa = R_Calloc(n_cell, double);
        double *tmpb = R_Calloc(n_cell, double);
        memset(tmpa, 0, sizeof(double)*n_cell);
        memset(tmpb, 0, sizeof(double)*n_cell);
        double mna = 0,
            mnb = 0;
        int j;
        for (j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            int cid = ai[j];
            tmpa[cid] = ax[j];
        }
        
        for (j = 0; j < n_cell; ++j) tmpb[j] = bx[j];

        for (j = 0; j < n_cell; ++j) {
            if (sens) {
                tmpb[j] = tmpb[j] - tmpa[j];
                if (tmpb[j] < 0) tmpb[j] = 0;
            }
            tmpb[j] = tmpb[j]/REAL(cs)[j];
            mnb += tmpb[j];
        }
        mnb = mnb/n_cell;

        for (j = ap[ii]; j < ap[ii+1]; ++j) {
            if (ISNAN(ax[j])) continue;            
            int cid = ai[j];
            tmpa[cid] = ax[j]/REAL(cs)[cid]*scale_factor;
            mna += tmpa[cid];
        }
        mna = mna/n_cell;
        
        double *tmpa_s = R_Calloc(n_cell, double);
        double *tmpb_s = R_Calloc(n_cell, double);
        smooth_W(tmpa, tmpa_s, n_cell, W);
        smooth_W(tmpb, tmpb_s, n_cell, W);
        double mna_s = 0,
            mnb_s = 0;
        for (j = 0; j < n_cell; ++j) {
            mna_s += tmpa_s[j];
            mnb_s += tmpb_s[j];
        }
        mna_s = mna_s/(double)n_cell;
        mnb_s = mnb_s/(double)n_cell;
        
        double Lx1 = 0,
            Lx2 = 0,
            //Ly1 = 0,
            //Ly2 = 0,
            ra  = 0,
            rb1 = 0,
            rb2 = 0;
        for (j = 0; j < n_cell; ++j) {
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
            shuffle(tmpa, ris[k], n_cell);
            smooth_W(tmpa, tmpa_s, n_cell, W);
            mna_s = 0;
            for (j = 0; j < n_cell; ++j) {
                mna_s += tmpa_s[j];
            }
            mna_s = mna_s/(double)n_cell;
            Lx1 = 0;
            ra = 0;
            rb1 = 0;
            for (j = 0; j < n_cell; ++j) {
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
*/
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
    // W = M_cholmod_transpose(W, (int)W->xtype, &c);
    
    int N = A->nrow; // cells
    int nc = A->ncol; // features

    double N1 = N -1;
    double N2 = N -2;
    double N3 = N -3;
    // double N4 = N -4;
    double NN = pow(N, 2);
    double tmp1 = NN - 3 *N + 3;
    double tmp2 = N1 *N2 * N3;
    double tmp3 = NN - N;

    double EI = -1/N1;
    double EI2 = pow(EI,2);
    
    SEXP Ival = PROTECT(allocVector(REALSXP, nc));
    //SEXP Cval = PROTECT(allocVector(REALSXP, nc));
    SEXP IZval = PROTECT(allocVector(REALSXP, nc));
    //SEXP CXval = PROTECT(allocVector(REALSXP, nc));    
    
    const int * ap = (int*)A->p;
    const int * ai = (int*)A->i;
    const double * ax = (double*)A->x;
    
    const int * bp = (int*)W->p;
    const int * bi = (int*)W->i;
    const double * bx = (double*)W->x;

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
    // const double varC_norm = ((2*S1+S2)*N1 - 4 * S02)/(2*(N+1)*S02);

    int ci;
#pragma omp parallel for num_threads(n_thread)
    for (ci = 0; ci < nc; ++ci) { // for each feature
        if (ap[ci] == ap[ci+1]) {
            REAL(Ival)[ci] = 0;
            continue;   
        }

        double *tmpa = R_Calloc(N, double);
        // memset(tmpa, 0, sizeof(double)*N);
        
        /* for (j = ap[ci]; j < ap[ci+1]; ++j) { */
        /*     if (ISNAN(ax[j])) continue; */
        /*     mn += ax[j]; */
        /* } */
        /* mn = mn/(double)N; */

        // scaled
        double mn = 0;
        for (int j = ap[ci]; j < ap[ci+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            tmpa[ai[j]] = ax[j];
            mn += ax[j];
        }
        mn = mn / N;

        double xix = 0;
        for (int j = 0; j < N; ++j) {
            xix += pow(tmpa[j] -mn, 2);
        }
        double sigma = sqrt(xix/N1);

        for (int j = 0; j < N; ++j) {
            tmpa[j] = (tmpa[j] -mn)/sigma;
        }

        // for (j = 0; j < N; ++j) tmpa[j] = tmpa[j] -mn;
        // xix = 0;
        double m2 = (double)N1/N;
        double m4 = 0;
        for (int j = 0; j < N; ++j) { // for cells
            //xix += pow(tmpa[j], 2);
            m4  += pow(tmpa[j], 4);
        }
        m4 = m4/N;
        
        // m2 = xix/(double)N;

        /* for (wi = 0; wi < N; ++wi) { */
        /*     if (bp[wi] == bp[wi+1]) continue; */
        /*     int wj; */
        /*     for (wj = bp[wi]; wj < bp[wi+1]; ++wj) { */
        /*         // Rprintf("W : %f, Zi : %f, Zj : %f\n", bx[wj], tmpa[wi], tmpa[bi[wj]]); */
        /*         C += bx[wj] * pow(tmpa[wi] - tmpa[bi[wj]],2); */
        /*     } */
        /* } */

        double varI = 0;
        if (rand) { // random
            //m2 = xix/N;
            //m4 = m4/N;
            double b2 = m4/pow(m2,2);
            varI = N*(tmp1*S1 - N*S2 + 3*S02) - b2*(tmp3*S1 - 2*N*S2 + 6*S02);
            varI = varI/(tmp2*S02) - EI2;
            //varC = N1*S1*(NN-3*N+3-N1*b2) - N1*S2*(NN+3*N-6-(NN-N+2)*b2)/4 + S02*(NN-3-N1*N1*b2);
            //varC = varC/(N*N2*N3*S02);
        }

        double xij = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = bp[i]; j < bp[i+1]; ++j) {
                if (ISNAN(bx[j])) continue;
                int cid = bi[j];
                if (i == cid) continue; // i != j
                //xij += (bx[j] * (tmpa[i]-mn) * (tmpa[cid]-mn));
                xij += bx[j] * tmpa[i] *tmpa[cid];
            }
        }

        double I = (N*xij)/(N1*S0);
        //Rprintf("N : %d, xij : %f, S0 : %f, mn : %f\n", N, xij, S0, mn);
        
        //Rprintf("wC : %f, sigma : %f, S0 : %f\n", C, sigma2, S0);
        // C = C/(2*sigma2*S0);
        //Rprintf("Cvar : %f, %f\n", varC, varC_norm);
        R_Free(tmpa);
#pragma omp critical
        {
            REAL(Ival)[ci] = I;
            // REAL(Cval)[ci] = C;
            REAL(IZval)[ci] = rand ? (I-EI)/sqrt(varI) : (I-EI)/sqrt(varI_norm);
            // REAL(CXval)[ci] = rand ? (C-1)/sqrt(varC) : (C-1)/sqrt(varC_norm);
        }
    }


    M_cholmod_free_sparse(&A, &c);
    // M_cholmod_free_sparse(&W, &c);
    
    SEXP ta = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ta, 0, Ival);
    // SET_VECTOR_ELT(ta, 1, Cval);
    SET_VECTOR_ELT(ta, 1, IZval);
    // SET_VECTOR_ELT(ta, 3, CXval);

    UNPROTECT(3);
    return ta;
}

SEXP moransi_perm_test(SEXP _A, SEXP _W, SEXP _scaled, SEXP _threads, SEXP _binarized, SEXP _perm)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP W = AS_CHM_SP__(_W);
    Rboolean scaled = asLogical(_scaled);
    int n_thread = asInteger(_threads);
    Rboolean binarized = asLogical(_binarized);
    int perm = asInteger(_perm);
    
    if (A->stype) return mkString("A cannot be symmetric");
    if (W->stype) return mkString("W cannot be symmetric");

    R_CheckStack();

    if (A->ncol != W->ncol) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix");
    if (A->ncol < 2) return mkString("Too few cells.");

    A = M_cholmod_transpose(A, (int)A->xtype, &c);
    //W = M_cholmod_transpose(W, (int)W->xtype, &c);
    
    int N = A->nrow;
    int nc = A->ncol; // features

    const int * ap = (int*)A->p;
    const int * ai = (int*)A->i;
    const double * ax = (double*)A->x;
    
    const int *bp = (int*)W->p;
    const int *bi = (int*)W->i;
    const double  *bx = (double*)W->x;

    double S0 = 0;
    
    int i;
    for (i = 0; i < W->nzmax; ++i) {
        if (bx[i] <= 0) continue;
        if (binarized) {
            S0++;
        } else {
            S0 += bx[i];
        }
    }

    SEXP Ival = PROTECT(allocVector(REALSXP, nc));
    SEXP Tval = PROTECT(allocVector(REALSXP, nc));

    int **ris = NULL;
    if (perm > 0) {
        ris = R_Calloc(perm, int*);
        for (int pi = 0; pi < perm; ++pi) {
            ris[pi] = random_idx(N);
        }
    }

    double t = 0;
    
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < nc; ++i) {
        if (ap[i] == ap[i+1]) {
            REAL(Ival)[i] = 0;
            REAL(Tval)[i] = 0;
            continue;
        }

        double *tmpa = R_Calloc(N, double);
        for (int j = ap[i]; j < ap[i+1]; ++j) {
            if (ISNAN(ax[j])) continue;                
            tmpa[ai[j]] = ax[j];
        }
            
        if (!scaled) {
            double mn = 0;
            double xix = 0;
            double sigma = 0;
            int N1 = N-1;
            int j;
            for (j = ap[i]; j < ap[i+1]; ++j) {
                if (ISNAN(ax[j])) continue;                
                mn += ax[j];
            }
            mn = mn / N;
            for (j = 0; j < N; ++j) {
                xix += pow(tmpa[j] -mn, 2);
            }
            sigma = sqrt(xix/N1);
            for (j = 0; j < N; ++j) {
                tmpa[j] = (tmpa[j] -mn)/sigma;
            }
        }
        
        double z2 = 0;
        for (int j = 0; j < N; ++j) {
            z2 += pow(tmpa[j],2);
        }

        double C = 0;
        for (int wi = 0; wi < N; ++wi) {
            if (bp[wi] == bp[wi+1]) continue;         
            for (int wj = bp[wi]; wj < bp[wi+1]; ++wj) {
                if (ISNAN(bx[wj])) continue;
                if (bx[wj] <= 0) continue;
                if (binarized) {
                    C += tmpa[wi] * tmpa[bi[wj]];
                } else {
                    C += bx[wj] * tmpa[wi] * tmpa[bi[wj]];
                }
            }
        }

        if (perm > 0) {
            double *Cs = R_Calloc(perm, double);
            double mn = 0;
            for (int k = 0; k < perm; ++k) {
                shuffle(tmpa, ris[k], N);
                Cs[k] = 0;
                for (int wi = 0; wi < N; ++wi) {
                    if (bp[wi] == bp[wi+1]) continue;         
                    for (int wj = bp[wi]; wj < bp[wi+1]; ++wj) {
                        if (ISNAN(bx[wj])) continue;
                        if (bx[wj] <= 0) continue;
                        if (binarized) {
                            Cs[k] += tmpa[wi] * tmpa[bi[wj]];
                        } else {
                            Cs[k] += bx[wj] * tmpa[wi] * tmpa[bi[wj]];
                        }
                    }
                }
                mn += Cs[k];
            }
            mn = mn/perm;
            double var = 0;
            for (int k = 0; k < perm; ++k) {
                var += pow((Cs[k]-mn),2);
            }
            var = sqrt(var/perm);
            
            t = (C - mn)/var;
            R_Free(Cs);
        }
        
        R_Free(tmpa);
     
#pragma omp critical
        {
            REAL(Ival)[i] = N/S0 * C/z2;
            REAL(Tval)[i] = t;
        }
    }

    if (perm > 0) {
        for (int pi = 0; pi < perm; ++pi) R_Free(ris[pi]);
        R_Free(ris);
    }

    M_cholmod_free_sparse(&A, &c);

    SEXP ta = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ta, 0, Ival);
    SET_VECTOR_ELT(ta, 1, Tval);

    UNPROTECT(3);
    return(ta);
}

