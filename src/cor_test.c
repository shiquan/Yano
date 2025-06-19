#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <assert.h>

#include "perm.h"

#ifndef Matrix_stubs
#define Matrix_stubs
#include "Matrix_stubs.c"
#endif

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

void smooth_W2(double * const a, double *s, const int N, CHM_SP W)
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
            s[i] += a[idx];
        }
    }
}

SEXP D_test_v1(SEXP _A,
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
    cholmod_common c;
    M_R_cholmod_start(&c);
            
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

    SEXP Dval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Tval  = PROTECT(allocVector(REALSXP, N_feature));
    
    const int * ap = (int*)A->p;
    const int * ai = (int*)A->i;
    const double * ax = (double*)A->x;

    const int * bp = (int*)B->p;
    const int * bi = (int*)B->i;
    const double * bx = (double*)B->x;

    random_index_init(perm, n_cell);
    
    R_CheckUserInterrupt();

    int i;
    
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < N_feature; ++i) {
        int ii = INTEGER(idx)[i]  -1;
        int ij = INTEGER(bidx)[i] -1;
        if (ap[ii] == ap[ii+1] || bp[ij] == bp[ij+1]) {
            REAL(Dval)[i]  = NA_REAL;
            REAL(Tval)[i]  = NA_REAL;
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
        
        for (int k = 0; k < perm; ++k) {
            if (method == 1 || method == 3) {
                shuffle_index(tmpa, k, n_cell);
                // shuffle(tmpa, ris[k], n_cell);
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
                shuffle_index(tmpb, k, n_cell);
                // shuffle(tmpb, ris[k], n_cell);
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
            REAL(Dval)[i]  = e;
            REAL(Tval)[i]  = t;
        }

        if (debug) {
            Rprintf("t, %f, mean, %f, var, %f \n", t, mean, var);
        }

    }

    random_index_free();
    
    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&B, &c);

    M_cholmod_finish(&c);
    SEXP ta = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ta, 0, Dval);
    SET_VECTOR_ELT(ta, 1, Tval);
    UNPROTECT(3);
    return ta;
}

extern CHM_SP imputation0(CHM_SP x, CHM_SP W, double filter, cholmod_common *c);

CHM_SP init_perm_matrix(CHM_SP X, int perm, cholmod_common *c)
{
    // M_R_cholmod_start(&c);    
    int *xi = (int*)X->i;
    int *xp = (int*)X->p;
    int i, p;
    for (i = 0; i < perm; ++i) { // start from 1, because 0 is orginal X
        for (p = xp[i+1]; p < xp[i+2]; ++p) {
            xi[p] = get_perm_idx(i, xi[p]);
        }
    }
    
    M_cholmod_sort(X, c);
    return X;
}

SEXP D_test_v2(SEXP _A,
               SEXP _B,
               SEXP _W,
               SEXP _permut,
               SEXP _threads,
               SEXP idx,
               SEXP bidx,
               SEXP _filter,
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

    double filter = asReal(_filter);
    const int perm = asInteger(_permut);
    int n_thread = asInteger(_threads);

    Rboolean debug = asLogical(_debug);    
    const int seed = asInteger(_seed);
    srand(seed);

    cholmod_common c;
    M_R_cholmod_start(&c);

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

    const int * ap = (int*)A->p;
    const int * ai = (int*)A->i;
    const double * ax = (double*)A->x;

    const int * bp = (int*)B->p;
    const int * bi = (int*)B->i;
    const double * bx = (double*)B->x;

    const int * wp = (int*)W->p;
    const int * wi = (int*)W->i;
    const double * wx = (double*)W->x;

    random_index_init(perm, n_cell);
    
    R_CheckUserInterrupt();

    SEXP Dval  = PROTECT(allocVector(REALSXP, N_feature));
    SEXP Tval  = PROTECT(allocVector(REALSXP, N_feature));
    
    int i;    
#pragma omp parallel for num_threads(n_thread)
    for (i = 0; i < N_feature; ++i) {
        cholmod_common c1;
        M_R_cholmod_start(&c1);
        int ii = INTEGER(idx)[i]  -1;
        int ij = INTEGER(bidx)[i] -1;
        if (ap[ii] == ap[ii+1] || bp[ij] == bp[ij+1]) {
            REAL(Dval)[i]  = NA_REAL;
            REAL(Tval)[i]  = NA_REAL;
            continue;
        }

        // init data for X and permutated Xs
        int fl = ap[ii+1] - ap[ii];
        int xnz = fl*(perm+1);
        CHM_SP XX = M_cholmod_allocate_sparse(n_cell, perm+1, xnz, FALSE, TRUE, 0, CHOLMOD_REAL, &c1);
        int j, s, p, k;
        int *xxp = (int*)XX->p;
        int *xxi = (int*)XX->i;
        double *xxx = (double*)XX->x;
        int n = 0;
        for (j = 0; j < perm+1; ++j) {
            xxp[j] = n;
            for (p = ap[ii]; p <ap[ii+1]; ++p, ++n) {
                xxi[n] = ai[p];
                xxx[n] = ax[p];
            }
        }
        xxp[j] = n;

        double mx = 0; // mean of X
        for (j = 0; j < xxp[1]; ++j) {
            mx += xxx[j];
        }
        mx = mx/n_cell;
        double X2 = 0;
        for (j = 0; j < xxp[1]; ++j) {
            X2 += pow(xxx[j] - mx,2);
        }

        XX = init_perm_matrix(XX, perm, &c1);
        CHM_SP tmp = XX;
        XX = M_cholmod_transpose(XX, (int)XX->xtype, &c1);
        M_cholmod_free_sparse(&tmp, &c1);
        
        CHM_SP SX = imputation0(XX, W, filter, &c1);
        M_cholmod_free_sparse(&XX, &c1);
        
        CHM_SP SXt = M_cholmod_transpose(SX, (int)SX->xtype, &c1);
        M_cholmod_free_sparse(&SX, &c1);
        
        double *tx = (double*)SXt->x;
        int *ti = (int*)SXt->i;
        int *tp = (int*)SXt->p;

        // smooth Y
        double Y[n_cell];
        memset(Y, 0, sizeof(double)*n_cell);
        double msy = 0; // mean value of smoothed Y (Y~)
        for (s = 0; s < n_cell; ++s) {
            p = bp[ij];
            for (j = wp[s]; j < wp[s+1]; ++j) {
                int id = wi[j];
                for (; p < bp[ij+1]; p++) {
                    if (bi[p] == id) {
                        Y[s] += bx[p]*wx[j];
                    } else if (bi[p] > id) {
                        break;
                    }
                }
                if (p == bp[ij+1]) break;
            }
            msy += Y[s];
        }
        
        // calculate Sx and r
        msy = msy/n_cell;
        double Yi[n_cell]; // Y~i - mean(Y~)
        double Yi2 = 0; // sum(Y~u - mean(Y~))^2
        for (s = 0; s < n_cell; ++s) {
            Yi[s] = Y[s] - msy;
            Yi2 += pow(Yi[s],2); 
        }
        Yi2 = sqrt(Yi2);

        if (debug) {
            Rprintf("X_mean, %f, smooth_Y_mean, %f, Yi^2, %f \n", mx, msy, Yi2);
        }

        double D[perm+1];
        double Xi[n_cell];
        for (j = 0; j < perm+1; ++j) {
            double Xi2 = 0;//, Xi3 = 0;
            double msx = 0;
            double SS = 0;
            memset(Xi, 0, sizeof(double)*n_cell);
            for (p = tp[j]; p < tp[j+1]; ++p) {
                msx += tx[p];
                Xi[ti[p]] = tx[p];
            }
            msx = msx/n_cell;
            double r = 0;
            for (p = 0; p < n_cell; ++p) {
                double xi0 = Xi[p] - msx;
                r += xi0*Yi[p];
                Xi2 += pow(xi0, 2);
            }
            SS = Xi2/X2;
            SS = sqrt(SS);
            Xi2 = sqrt(Xi2);
            r = r / (Xi2 * Yi2);
            D[j] = SS * (1-r);

            if (debug && j == 0) {
                Rprintf("smooth_X_mean, %f, Xi^2, %f, SS, %f, r, %f, D, %f \n", msx, Xi2, SS, r, D[j]);
            }
        }

        M_cholmod_free_sparse(&SXt, &c1);
        M_cholmod_finish(&c1);
        // calculate mean, var
        double md = 0, var = 0;        
        for (j = 1; j < perm+1; ++j) {
            md += D[j];
        }
        md = md/perm;
        for (j = 1; j < perm+1; ++j) {
            var += pow(D[j]-md, 2);
        }
        var = sqrt(var/perm);

        double t = (D[0] - md)/var;
        
/* #pragma omp critical */
/*         { */
        REAL(Dval)[i]  = D[0];
        REAL(Tval)[i]  = t;
        /* } */

        if (debug) {
            Rprintf("t, %f, mean, %f, var, %f \n", t, md, var);
        }
    } // end of loop, for each feature

    random_index_free();

    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&B, &c);
    M_cholmod_finish(&c);
    R_CheckStack();
    
    SEXP ta = PROTECT(allocVector(VECSXP, 2));

    SET_VECTOR_ELT(ta, 0, Dval);
    SET_VECTOR_ELT(ta, 1, Tval);

    UNPROTECT(3);
    return ta;
}

SEXP D_distribution_test_v2(SEXP _A,
                            SEXP _B,
                            SEXP _W,
                            SEXP _permut,
                            SEXP _filter,
                            SEXP _seed,
                            SEXP _debug)
{
    int la = length(_A);
    int lb = length(_B);
    
    CHM_SP W = AS_CHM_SP__(_W);
    double filter = asReal(_filter);
    const int perm = asInteger(_permut);
    
    if (W->stype) return mkString("W cannot be symmetric");
    if (la != lb) return mkString("Unequal length of A and B");
   
    const int n_cell = W->nrow;
    
    if (la != n_cell) return mkString("Unequal length of A and W");
    
    Rboolean debug = asLogical(_debug);    
    const int seed = asInteger(_seed);
    srand(seed);

    cholmod_common c;
    M_R_cholmod_start(&c);
    
    const int *wp = (int*)W->p;
    const int *wi = (int*)W->i;
    const double *wx = (double*)W->x;

    random_index_init(perm, n_cell);

    int i, j, s, p, k;
    double X[n_cell];
    double Y[n_cell];
    double mx = 0;
    int nl = 0;
    for (i = 0; i < n_cell; ++i) {
        X[i] = REAL(_A)[i];
        Y[i] = REAL(_B)[i];
        if (X[i] == 0) continue;
        nl++;
        mx += X[i];
    }
    mx = mx/n_cell;
    Rprintf("mx, %f\n", mx);
    double X2 = 0;
    for (i = 0; i < n_cell; ++i) {
        X2 += pow(X[i] - mx,2);
    }
    Rprintf("X2, %f\n", X2);
    int xnz = nl*(perm+1);
    CHM_SP XX = M_cholmod_allocate_sparse(n_cell, perm+1,xnz, FALSE, TRUE, 0, CHOLMOD_REAL, &c);

    int *xxp = (int*)XX->p;
    int *xxi = (int*)XX->i;
    double *xxx = (double*)XX->x;
    int n = 0;
    for (i = 0; i < perm+1; ++i) {
        xxp[i] = n;
        for (j = 0; j < n_cell; ++j) {
            if (X[j] == 0) continue;
            xxi[n] = j;
            xxx[n] = X[j];
            n++;
        }
    }

    xxp[i] = n;
    assert(n == xnz);

    Rprintf("init perm\n");
    XX = init_perm_matrix(XX, perm, &c);

    CHM_SP tmp = XX;
    XX = M_cholmod_transpose(XX, (int)XX->xtype, &c);
    M_cholmod_free_sparse(&tmp, &c);

    CHM_SP SX = imputation0(XX, W, filter, &c);
    M_cholmod_free_sparse(&XX, &c);
    CHM_SP SXt = M_cholmod_transpose(SX, (int)SX->xtype, &c);
    M_cholmod_free_sparse(&SX, &c);

    double *tx = (double*)SXt->x;
    int *ti = (int*)SXt->i;
    int *tp = (int*)SXt->p;

    double Ys[n_cell];
    memset(Ys, 0, sizeof(double)*n_cell);
    double msy = 0; // mean value of smoothed Y (Y~)
    for (s = 0; s < n_cell; ++s) {
        for (j = wp[s]; j < wp[s+1]; ++j) {
            int id = wi[j];
            if (Y[id] == 0) continue;
            Ys[s] += Y[id]*wx[j];
        }
        msy += Ys[s];
    }
    // calculate Sx and r
    msy = msy/n_cell;
    double Yi[n_cell]; // Y~i - mean(Y~)
    double Yi2 = 0; // sum(Y~u - mean(Y~))^2
    for (s = 0; s < n_cell; ++s) {
        Yi[s] = Ys[s] - msy;
        Yi2 += pow(Yi[s],2); 
    }
    Yi2 = sqrt(Yi2);

    if (debug) {
        Rprintf("X_mean, %f, smooth_Y_mean, %f, Yi^2, %f \n", mx, msy, Yi2);
    }

    double D[perm+1];
    double Xi[n_cell];
    for (j = 0; j < perm+1; ++j) {
        double Xi2 = 0;//, Xi3 = 0;
        double msx = 0;
        double SS = 0;
        memset(Xi, 0, sizeof(double)*n_cell);
        for (p = tp[j]; p < tp[j+1]; ++p) {
            msx += tx[p];
            Xi[ti[p]] = tx[p];
        }
        msx = msx/n_cell;
        double r = 0;
        for (p = 0; p < n_cell; ++p) {
            double xi0 = Xi[p] - msx;
            //double xi1 = Xi[p] - mx;
            r += xi0*Yi[p];
            Xi2 += pow(xi0, 2);
            //Xi3 += pow(xi1, 2);
        }
        SS = Xi2/X2;
        SS = sqrt(SS);
        Xi2 = sqrt(Xi2);
        r = r / (Xi2 * Yi2);            
        D[j] = SS * (1-r);
        
        if (debug && j == 0) {
            Rprintf("smooth_X_mean, %f, Xi^2, %f, SS, %f, r, %f, D, %f \n", msx, Xi2, SS, r, D[j]);
        }
    }
    
    M_cholmod_free_sparse(&SXt, &c);
    M_cholmod_finish(&c);
    SEXP DD = PROTECT(allocVector(REALSXP, perm+1));
    for (i = 0; i < perm+1; ++i) {
        REAL(DD)[i] = D[i];
    }

    UNPROTECT(1);
    return DD;
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

SEXP D_distribution_test2(SEXP _A, SEXP _B, SEXP _W, SEXP _permut, SEXP _threads)
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
    smooth_W2(b, tmpb_s, N, W);

    double mnb_s = 0;
    for (i = 0; i < N; ++i) {
        mnb_s += tmpb_s[i];
    }
    mnb_s = mnb_s/N;

    double *tmpa_s = R_Calloc(N, double);
    smooth_W2(a, tmpa_s, N, W);
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
        smooth_W2(s, tmpa_s, N, W);
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

