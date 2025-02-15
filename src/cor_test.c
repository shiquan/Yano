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
        double e;
        if (method == 1) {
            e = sqrt(Lx) * (1-r);
        } else if (method == 2) {
            // e = sqrt(Lx) * sqrt(Ly) * (1-r);
            e = sqrt(Ly) * (1-r);
        } else {
            e = sqrt(Lx) * sqrt(Ly) * r;
        }
        
        if (debug) {
            Rprintf("Lx, %f, D, %f, r, %f\n", Lx, e, r);
        }

        double r0 = r;
        double e0 = e;

        // mean, var
        double mean = 0, var = 0;
        double *es = R_Calloc(perm, double);
        // double *es = alloca(perm *sizeof(double));
        // R_CheckStack();
        
        for (int k = 0; k < perm; ++k) {
            if (method == 1 || method == 3) {
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
            }
            if (method == 2) {
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
                // es[k] = sqrt(Lx) * sqrt(Ly) *(1-r);
                es[k] = sqrt(Ly) *(1-r);
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
            REAL(Rval)[i]  = r;
            REAL(Dval)[i]  = e;
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
