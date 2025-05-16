#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>

extern int *random_idx(const int n);
extern void shuffle(double tmp[], int const idx[], const int n);

SEXP association_test1(SEXP A, SEXP B, SEXP C, SEXP _W, SEXP _perm, SEXP return_arr)
{
    int n_perm = asInteger(_perm);
    Rboolean ret = asLogical(return_arr);
    
    CHM_SP W = AS_CHM_SP__(_W);

    if (W->stype) return mkString("W cannot be symmetric");

    int la = length(A);
    int lb = length(B);
    int lc = length(C);
    if (la != lb || la != lc || lb != lc) {
        mkString("Unequal length of input array.");
    }
    if (la != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (W->ncol < 2) return mkString("Too few cells."); // to do
    
    const int * wp = (int*)W->p;
    const int * wi = (int*)W->i;
    const double * wx = (double*)W->x;

    R_CheckUserInterrupt();

    int n_cell = la;
    double tmpa[n_cell];
    double tmpb[n_cell];
    double tmpc[n_cell];
    
    memset(tmpa, 0, sizeof(double)*n_cell);
    memset(tmpb, 0, sizeof(double)*n_cell);
    memset(tmpc, 0, sizeof(double)*n_cell);

    double mna = 0, mnb = 0, mnc = 0;

    int i;
    int **ris = NULL;
    ris = R_Calloc(n_perm,int*);
    for (i = 0; i < n_perm; ++i) {
        ris[i] = random_idx(n_cell);
    }

    for (i = 0; i < n_cell; ++i) {
        tmpa[i] = REAL(A)[i];
        tmpb[i] = REAL(B)[i];
        tmpc[i] = REAL(C)[i];
        mna += REAL(A)[i];
        mnb += REAL(B)[i];
        mnc += REAL(C)[i];
    }
    mna = mna/n_cell;
    mnb = mnb/n_cell;
    mnc = mnc/n_cell;
    
    double sd1 = 0;
    double sd2 = 0;
    double sd3 = 0;
    int j;
    for (j = 0; j < n_cell; ++j) {
        tmpa[j] = tmpa[j] - mna;
        tmpb[j] = tmpb[j] - mnb;
        tmpc[j] = tmpc[j] - mnc;
        sd1 += tmpa[j] * tmpa[j];
        sd2 += tmpb[j] * tmpb[j];
        sd3 += tmpc[j] * tmpc[j];
    }
    sd1 = sqrt(sd1/n_cell);
    sd2 = sqrt(sd2/n_cell);
    sd3 = sqrt(sd3/n_cell);

    double I1 = 0.0;
    double I2 = 0.0;
    int ii;
    int jj;
    for (ii = 0; ii < n_cell; ++ii) {
        if (tmpa[ii] == 0) continue;
        for (jj = wp[ii]; jj < wp[ii+1]; ++jj) {
            if (ISNAN(wx[jj])) continue;
            if (wx[jj] <= 0) continue;
            int idx = wi[jj];
            I1 += wx[jj] * tmpa[ii] * tmpb[idx];
            I2 += wx[jj] * tmpa[ii] * tmpc[idx];
        }
    }

    I1 = I1/(n_cell * sd1 * sd2);
    I2 = I2/(n_cell * sd1 * sd3);

    double delta = I1 - I2;
    SEXP ds = PROTECT(allocVector(REALSXP, n_perm));

    for (i = 0; i < n_perm; ++i) {
        shuffle(tmpa, ris[i], n_cell);
        int ii, jj;
        double i1 = 0, i2 = 0;
        for (ii = 0; ii < n_cell; ++ii) {
            if (tmpa[ii] == 0) continue;
            for (jj = wp[ii]; jj < wp[ii+1]; ++jj) {
                if (ISNAN(wx[jj])) continue;
                if (wx[jj] <= 0) continue;
                int idx = wi[jj];
                i1 += wx[jj] * tmpa[ii] * tmpb[idx];
                i2 += wx[jj] * tmpa[ii] * tmpc[idx];
            }
        }
        i1 = i1/(n_cell * sd1 * sd2);
        i2 = i2/(n_cell * sd1 * sd3);
        REAL(ds)[i] = i1 - i2;
    }
    
    for (j = 0; j < n_perm; ++j) R_Free(ris[j]);
    if (n_perm) R_Free(ris);

    if (ret) {
        SEXP ta = PROTECT(allocVector(VECSXP, 3));
        SET_VECTOR_ELT(ta, 0, ScalarReal(I1));
        SET_VECTOR_ELT(ta, 1, ScalarReal(I2));
        SET_VECTOR_ELT(ta, 2, ds);
        UNPROTECT(2);
        return ta;
    } else {
        double mn = 0, sd = 0;
        for (i =0; i < n_perm; ++i) {
            mn += REAL(ds)[i];
        }
        for (i = 0; i < n_perm; ++i) {
            REAL(ds)[i] = REAL(ds)[i] - mn;
            sd += REAL(ds)[i] *REAL(ds)[i];
        }
        sd = sqrt(sd/n_perm);
        
        SEXP ta = PROTECT(allocVector(VECSXP, 3));
        SET_VECTOR_ELT(ta, 0, ScalarReal(I1));
        SET_VECTOR_ELT(ta, 1, ScalarReal(I2));
        SET_VECTOR_ELT(ta, 2, ScalarReal((I1-I2-mn)/sd));
        UNPROTECT(2);
        return ta;
    }
}

SEXP association_test2(SEXP _A, SEXP B, SEXP C, SEXP _W, SEXP _perm, SEXP _threads)
{
    int n_thread = asInteger(_threads);
    int n_perm = asInteger(_perm);

    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP W = AS_CHM_SP__(_W);

    if (W->stype) return mkString("W cannot be symmetric.");

    int la = A->nrow;
    int lb = length(B);
    int lc = length(C);
    if (la != lb || la != lc || lb != lc) {
        mkString("Unequal length of input array.");
    }
    if (la != W->nrow) return mkString("A column and W row do not match.");
    if (W->nrow != W->ncol) return mkString("W is not a square matrix.");
    if (W->ncol < 2) return mkString("Too few cells."); // to do
    
    const int *wp = (int*)W->p;
    const int *wi = (int*)W->i;
    const double *wx = (double*)W->x;

    const int *ap = (int*)A->p;
    const int *ai = (int*)A->i;
    const double *ax = (double*)A->x;

    int n_feature = A->ncol;
    if (n_thread > n_feature) n_thread = n_feature;

    int n_cell = la;

    SEXP I1  = PROTECT(allocVector(REALSXP, n_feature));
    SEXP I2  = PROTECT(allocVector(REALSXP, n_feature));
    SEXP T  = PROTECT(allocVector(REALSXP, n_feature));
        
    double tmpb[n_cell];
    double tmpc[n_cell];
    memset(tmpb, 0, sizeof(double)*n_cell);
    memset(tmpc, 0, sizeof(double)*n_cell);

    double mnb = 0, mnc = 0;

    int i;
    int **ris = NULL;
    ris = R_Calloc(n_perm,int*);
    for (i = 0; i < n_perm; ++i) {
        ris[i] = random_idx(n_cell);
    }

    for (i = 0; i < n_cell; ++i) {
        tmpb[i] = REAL(B)[i];
        tmpc[i] = REAL(C)[i];
        mnb += REAL(B)[i];
        mnc += REAL(C)[i];
    }
    mnb = mnb/n_cell;
    mnc = mnc/n_cell;
    
    double sd2 = 0;
    double sd3 = 0;
    int j;
    for (j = 0; j < n_cell; ++j) {
        tmpb[j] = tmpb[j] - mnb;
        tmpc[j] = tmpc[j] - mnc;
        sd2 += tmpb[j] * tmpb[j];
        sd3 += tmpc[j] * tmpc[j];
    }
    sd2 = sqrt(sd2/n_cell);
    sd3 = sqrt(sd3/n_cell);

    R_CheckUserInterrupt();

    int k;
#pragma omp parallel for num_threads(n_thread)    
    for (k = 0; k < n_feature; ++k) {
        double tmpa[n_cell];
        memset(tmpa, 0, sizeof(double)*n_cell);
        double sd1 = 0;
        int i, j;
        double mna = 0;
        for (i = ap[k]; i < ap[k+1]; ++i) {
            int idx = ai[i];
            tmpa[idx] = ax[i];
            mna += ax[i];
        }
        mna = mna/n_cell;
        for (j = 0; j < n_cell; ++j) {
            tmpa[j] = tmpa[j] - mna;
            sd1 += tmpa[j] * tmpa[j];
        }
        sd1 = sqrt(sd1/n_cell);

        double i1 = 0;
        double i2 = 0;
        int ii;
        int jj;
        for (ii = 0; ii < n_cell; ++ii) {
            if (tmpa[ii] == 0) continue;
            for (jj = wp[ii]; jj < wp[ii+1]; ++jj) {
                if (ISNAN(wx[jj])) continue;
                if (wx[jj] <= 0) continue;
                int idx = wi[jj];
                i1 += wx[jj] * tmpa[ii] * tmpb[idx];
                i2 += wx[jj] * tmpa[ii] * tmpc[idx];
            }
        }
        
        i1 = i1/(n_cell * sd1 * sd2);
        i2 = i2/(n_cell * sd1 * sd3);
        
        double delta = i1 - i2;
        double ds[n_perm];
        memset(ds, 0, sizeof(double)*n_perm);
        double mn = 0, sd = 0;
        for (i = 0; i < n_perm; ++i) {
            shuffle(tmpa, ris[i], n_cell);
            int ii, jj;
            double i1 = 0, i2 = 0;
            for (ii = 0; ii < n_cell; ++ii) {
                if (tmpa[ii] == 0) continue;
                for (jj = wp[ii]; jj < wp[ii+1]; ++jj) {
                    if (ISNAN(wx[jj])) continue;
                    if (wx[jj] <= 0) continue;
                    int idx = wi[jj];
                    i1 += wx[jj] * tmpa[ii] * tmpb[idx];
                    i2 += wx[jj] * tmpa[ii] * tmpc[idx];
                }
            }
            i1 = i1/(n_cell * sd1 * sd2);
            i2 = i2/(n_cell * sd1 * sd3);
            ds[i] = i1 - i2;
            mn += ds[i];
        }
        mn = mn/n_cell;

        for (i = 0; i < n_perm; ++i) {
            ds[i] = ds[i] - mn;
            sd += ds[i] * ds[i];
        }
        sd = sqrt(sd/n_perm);

        REAL(I1)[k] = i1;
        REAL(I2)[k] = i2;
        REAL(T)[k] = (delta-mn)/sd;
        /* REAL(M)[k] = mn; */
        /* REAL(S)[k] = sd; */
    }

    for (j = 0; j < n_perm; ++j) R_Free(ris[j]);
    if (n_perm) R_Free(ris);

    SEXP ta = PROTECT(allocVector(VECSXP, 3));
    SET_VECTOR_ELT(ta, 0, I1);
    SET_VECTOR_ELT(ta, 1, I2);
    SET_VECTOR_ELT(ta, 2, T);
    /* SET_VECTOR_ELT(ta, 3, M); */
    /* SET_VECTOR_ELT(ta, 4, S); */
    UNPROTECT(4);
    return ta;
}
