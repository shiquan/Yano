#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>
#include <assert.h>



#define CACHE_PER_BATCH 10000
#define MIN_HIT 1

extern int *random_idx(const int n);
extern void shuffle(double tmp[], int const idx[], const int n);

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
    cholmod_common c;
    M_R_cholmod_start(&c);
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

    int ci;
#pragma omp parallel for num_threads(n_thread)
    for (ci = 0; ci < nc; ++ci) { // for each feature
        if (ap[ci] == ap[ci+1]) {
            REAL(Ival)[ci] = 0;
            continue;   
        }

        double *tmpa = R_Calloc(N, double);
        
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

        double m2 = (double)N1/N;
        double m4 = 0;
        for (int j = 0; j < N; ++j) { // for cells
            m4  += pow(tmpa[j], 4);
        }
        m4 = m4/N;
        
        double varI = 0;
        if (rand) { // random
            double b2 = m4/pow(m2,2);
            varI = N*(tmp1*S1 - N*S2 + 3*S02) - b2*(tmp3*S1 - 2*N*S2 + 6*S02);
            varI = varI/(tmp2*S02) - EI2;
        }

        double xij = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = bp[i]; j < bp[i+1]; ++j) {
                if (ISNAN(bx[j])) continue;
                int cid = bi[j];
                if (i == cid) continue; // i != j
                xij += bx[j] * tmpa[i] *tmpa[cid];
            }
        }

        double I = (N*xij)/(N1*S0);
        R_Free(tmpa);
#pragma omp critical
        {
            REAL(Ival)[ci] = I;
            REAL(IZval)[ci] = rand ? (I-EI)/sqrt(varI) : (I-EI)/sqrt(varI_norm);
        }
    }

    M_cholmod_free_sparse(&A, &c);
    M_cholmod_finish(&c);
    SEXP ta = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ta, 0, Ival);
    SET_VECTOR_ELT(ta, 1, IZval);

    UNPROTECT(3);
    return ta;
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
    
    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (A->ncol < 2) return mkString("Too few cells."); // to do

    cholmod_common c;
    M_R_cholmod_start(&c);

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

    M_cholmod_finish(&c);

    SEXP ta = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(ta, 0, Ival);
    SET_VECTOR_ELT(ta, 1, Pval);

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

    cholmod_common c;
    M_R_cholmod_start(&c);        

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

    M_cholmod_finish(&c);

    SEXP ta = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ta, 0, Ival);
    SET_VECTOR_ELT(ta, 1, Tval);

    UNPROTECT(3);
    return(ta);
}

