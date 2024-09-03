#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>

extern int *random_idx(const int n);
extern void shuffle(double tmp[], int const idx[], const int n);
extern void smooth_W(double * const a, double *s, const int N, CHM_SP W);

struct val {
    int v1; // feature 
    int v2; // test.feature
    double L;
    double t;
    double r;
    double d;
    double t_d;
};
SEXP dl_score(SEXP _A, SEXP _B, SEXP _W, SEXP _perm, SEXP _threads, SEXP _seed)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP B = AS_CHM_SP__(_B);
    CHM_SP W = AS_CHM_SP__(_W);

    const int perm = asInteger(_perm);
    int n_thread = asInteger(_threads);

    if (A->nrow != B->nrow) return mkString("Unequal row of A and B.");
    if (A->nrow != W->nrow) return mkString("A row and W row do not match.");
    if (W->ncol != W->nrow) return mkString("W is not a square matrix.");

    const int *ap = (int *)A->p;
    const int *ai = (int *)A->i;
    const double *ax = (double *)A->x;

    const int *bp = (int *)B->p;
    const int *bi = (int *)B->i;
    const double *bx = (double *)B->x;

    const int n_cell = A->nrow;

    R_CheckStack();
    R_CheckUserInterrupt();

    const int seed = asInteger(_seed);
    srand(seed);

    int **ris = NULL;
    ris = R_Calloc(perm, int*);
    for (int pi = 0; pi < perm; ++pi) {
        ris[pi] = random_idx(n_cell);
    }

    int nn = 0;
    
    const int n_feature1 = A->ncol;
    const int n_feature2 = B->ncol;

    int n_record = n_feature1 * n_feature2;
    
    SEXP val1 = PROTECT(allocVector(INTSXP, n_record)); // feature
    SEXP val2 = PROTECT(allocVector(INTSXP, n_record)); // test.feature
    SEXP val3 = PROTECT(allocVector(REALSXP, n_record)); // L
    SEXP val4 = PROTECT(allocVector(REALSXP, n_record)); // t
    SEXP val5 = PROTECT(allocVector(REALSXP, n_record)); // Pearson's r
    SEXP val6 = PROTECT(allocVector(REALSXP, n_record)); // d
    SEXP val7 = PROTECT(allocVector(REALSXP, n_record)); // t_d
    
    if (n_feature1 >= n_feature2) {
        void **dat = R_Calloc(n_feature1, void*);
        
#pragma omp parallel for num_threads(n_thread)
        for (int i = 0; i < n_feature1; ++i) {
            struct val *vs = R_Calloc(n_feature2, struct val);
            for (int j = 0; j < n_feature2; ++j) {

                struct val *v = &vs[j];
                v->v1 = i;
                v->v2 = j;
                
                double mna = 0, mnb = 0;
                double *tmpa = R_Calloc(n_cell, double);
                double *tmpb = R_Calloc(n_cell, double);
                memset(tmpa, 0, sizeof(double)*n_cell);
                memset(tmpb, 0, sizeof(double)*n_cell);
                
                for (int ii = ap[i]; ii < ap[i+1]; ++ii) {
                    if (ISNAN(ax[ii])) continue;
                    int cid = ai[ii];
                    tmpa[cid] = ax[ii];
                    mna += tmpa[cid];
                }
                
                for (int ii = bp[j]; ii < bp[j+1]; ++ii) {
                    if (ISNAN(bx[ii])) continue;
                    int cid = bi[ii];
                    tmpb[cid] = bx[ii];
                    mnb += tmpb[cid];
                }
                
                mna = mna/n_cell;
                mnb = mnb/n_cell;

                double *tmpa_s = R_Calloc(n_cell, double);
                double *tmpb_s = R_Calloc(n_cell, double);
                smooth_W(tmpa, tmpa_s, n_cell, W);
                smooth_W(tmpb, tmpb_s, n_cell, W);
                double mna_s = 0;
                double mnb_s = 0;
                for (int k = 0; k < n_cell; ++k) {
                    mna_s += tmpa_s[k];
                    mnb_s += tmpb_s[k];
                }
                mna_s = mna_s/(double)n_cell;
                mnb_s = mnb_s/(double)n_cell;

                double Lx1 = 0, Lx2 = 0, Ly1 = 0, Ly2 = 0, ra = 0, rb1 = 0, rb2 = 0;
                double pr1 = 0;
                for (int k = 0; k < n_cell; ++k) {
                    Lx1 += pow(tmpa_s[k]-mna,2);
                    Lx2 += pow(tmpa[k]-mna,2);
                    Ly1 += pow(tmpb_s[k]-mnb,2);
                    Ly2 += pow(tmpb[k]-mnb,2);

                    pr1 += (tmpa[k]-mna)*(tmpb[k]-mnb);
                    tmpa_s[k] = tmpa_s[k] - mna_s;
                    tmpb_s[k] = tmpb_s[k] - mnb_s;
                    ra  += tmpa_s[k] * tmpb_s[k];
                    ra  += tmpa_s[k] * tmpb[k];
                    rb1 += pow(tmpa_s[k],2);
                    rb2 += pow(tmpb_s[k],2);
                }

                double pr = pr1/sqrt(Lx2*Ly2);
                rb1 = sqrt(rb1);
                rb2 = sqrt(rb2);
                
                double Lx = Lx1/Lx2;
                double Ly = Ly1/Ly2;
                double r = ra/(rb1*rb2);
                double L = sqrt(Lx) * sqrt(Ly) *r ;
                double D = sqrt(Lx) * (1-r);
                // Rprintf("mna = %f, mnb = %f, mna_s = %f, mnb_s = %f, Lx = %f, Ly = %f, r = %f, L = %f\n", mna, mnb, mna_s, mnb_s, Lx, Ly, r, L);
                        
                double mn = 0, var = 0;
                double mn_d = 0, var_d = 0;
                double *Ls = R_Calloc(perm, double);
                double *Ds = R_Calloc(perm, double);
                
                for (int k = 0; k < perm; ++k) {
                    shuffle(tmpa, ris[k], n_cell);
                    shuffle(tmpb, ris[k], n_cell);

                    smooth_W(tmpa, tmpa_s, n_cell, W);
                    smooth_W(tmpb, tmpb_s, n_cell, W);
                    
                    mna_s = 0;
                    mnb_s = 0;
                    for (int ii = 0; ii < n_cell; ++ii) {
                        mna_s += tmpa_s[ii];
                        mnb_s += tmpb_s[ii];
                    }
                    mna_s = mna_s/n_cell;
                    mnb_s = mnb_s/n_cell;
                    
                    Lx1 = Ly1 = ra = rb1 = rb2 = 0;                
                    for (int jj = 0; jj < n_cell; ++jj) {
                        Lx1 += pow(tmpa_s[jj] - mna, 2);
                        Ly1 += pow(tmpb_s[jj] - mnb, 2);
                        
                        tmpa_s[jj] = tmpa_s[jj] - mna_s;
                        tmpb_s[jj] = tmpb_s[jj] - mnb_s;
                        
                        ra += tmpa_s[jj] * tmpb_s[jj];
                        rb1 += pow(tmpa_s[jj],2);
                        rb2 += pow(tmpb_s[jj],2);
                    }
                    
                    rb1 = sqrt(rb1);
                    rb2 = sqrt(rb2);
                    r = ra/(rb1*rb2);
                    
                    Lx = Lx1/Lx2;
                    Ly = Ly1/Ly2;
                    
                    Ls[k] = sqrt(Lx) * sqrt(Ly) * r;
                    Ds[k] = sqrt(Lx) * (1-r);
                    mn += Ls[k];
                    mn_d += Ds[k];
                }
                mn = mn/perm;
                mn_d = mn_d/perm;
                
                for (int jj = 0; jj < perm; ++jj) {
                    var += pow((Ls[jj]-mn), 2);
                    var_d += pow((Ds[jj]-mn_d), 2);
                }
                var = sqrt(var/perm);
                var_d = sqrt(var_d/perm);
                
                double t = (L-mn)/var;
                double t_d = (D-mn_d)/var_d;
                v->t = t;
                v->L = L;
                v->r = pr;
                v->d = D;
                v->t_d = t_d;
                
                R_Free(Ls);
                R_Free(tmpa_s);
                R_Free(tmpb_s);
                R_Free(tmpa);
                R_Free(tmpb);
                R_Free(Ds);
            }
            
#pragma omp critical
            {
                dat[i] = vs;
            }
        }

        //
        for (int i = 0; i < n_feature1; ++i) {
            struct val *vs = dat[i];
            for (int j = 0; j < n_feature2; ++j) {
                struct val *v = &vs[j];
                INTEGER(val1)[nn] = v->v1+1;
                INTEGER(val2)[nn] = v->v2+1;
                REAL(val3)[nn] = v->L;
                REAL(val4)[nn] = v->t;
                REAL(val5)[nn] = v->r;
                REAL(val6)[nn] = v->d;
                REAL(val7)[nn] = v->t_d;                
                nn++;
            }
            R_Free(vs);
        }
        R_Free(dat);
    } else {
        void **dat = R_Calloc(n_feature2, void*);
                
#pragma omp parallel for num_threads(n_thread)
        for (int i = 0; i < n_feature2; ++i) {
            struct val *vs = R_Calloc(n_feature1, struct val);
            for (int j = 0; j < n_feature1; ++j) {
                struct val *v = &vs[j];
                v->v1 = j;
                v->v2 = i;
                
                double mna = 0, mnb = 0;
                double *tmpa = R_Calloc(n_cell, double);
                double *tmpb = R_Calloc(n_cell, double);
                memset(tmpa, 0, sizeof(double)*n_cell);
                memset(tmpb, 0, sizeof(double)*n_cell);
                
                for (int ii = ap[j]; ii < ap[j+1]; ++ii) {
                    if (ISNAN(ax[ii])) continue;
                    int cid = ai[ii];
                    tmpa[cid] = ax[ii];
                    mna += tmpa[cid];
                }
                
                for (int ii = bp[i]; ii < bp[i+1]; ++ii) {
                    if (ISNAN(bx[ii])) continue;
                    int cid = bi[ii];
                    tmpb[cid] = bx[ii];
                    mnb += tmpb[cid];
                }
                
                mna = mna/n_cell;
                mnb = mnb/n_cell;
                
                double *tmpa_s = R_Calloc(n_cell, double);
                double *tmpb_s = R_Calloc(n_cell, double);
                smooth_W(tmpa, tmpa_s, n_cell, W);
                smooth_W(tmpb, tmpb_s, n_cell, W);
                double mna_s = 0;
                double mnb_s = 0;
                for (int k = 0; k < n_cell; ++k) {
                    mna_s += tmpa_s[k];
                    mnb_s += tmpb_s[k];
                }
                mna_s = mna_s/n_cell;
                mnb_s = mnb_s/n_cell;
                
                double Lx1 = 0, Lx2 = 0, Ly1 = 0, Ly2 = 0, ra = 0, rb1 = 0, rb2 = 0;
                double pr1 = 0;
                for (int k = 0; k < n_cell; ++k) {
                    Lx1 += pow(tmpa_s[k]-mna,2);
                    Lx2 += pow(tmpa[k]-mna,2);
                    Ly1 += pow(tmpb_s[k]-mnb,2);
                    Ly2 += pow(tmpb[k]-mnb,2);

                    pr1 += (tmpa[k]-mna)*(tmpb[k]-mnb);
                    tmpa_s[k] = tmpa_s[k] - mna_s;
                    tmpb_s[k] = tmpb_s[k] - mnb_s;
                    ra  += tmpa_s[k] * tmpb_s[k];
                    ra  += tmpa_s[k] * tmpb[k];
                    rb1 += pow(tmpa_s[k],2);
                    rb2 += pow(tmpb_s[k],2);
                }
                
                rb1 = sqrt(rb1);
                rb2 = sqrt(rb2);

                double pr = pr1/sqrt(Lx2*Ly2);
                double Lx = Lx1/Lx2;
                double Ly = Ly1/Ly2;
                double r = ra/(rb1*rb2);
                double L = sqrt(Lx) * sqrt(Ly) *r ;
                double D = sqrt(Lx) * (1-r);
                // Rprintf("mna = %f, mnb = %f, mna_s = %f, mnb_s = %f, Lx = %f, Ly = %f, r = %f, L = %f\n", mna, mnb, mna_s, mnb_s, Lx, Ly, r, L);

                double mn = 0, var = 0;
                double mn_d = 0, var_d = 0;
                double *Ls = R_Calloc(perm, double);
                double *Ds = R_Calloc(perm, double);
                
                for (int k = 0; k < perm; ++k) {
                    shuffle(tmpa, ris[k], n_cell);
                    shuffle(tmpb, ris[k], n_cell);

                    smooth_W(tmpa, tmpa_s, n_cell, W);
                    smooth_W(tmpb, tmpb_s, n_cell, W);
                    
                    mna_s = 0;
                    mnb_s = 0;
                    for (int ii = 0; ii < n_cell; ++ii) {
                        mna_s += tmpa_s[ii];
                        mnb_s += tmpb_s[ii];
                    }
                    mna_s = mna_s/(double)n_cell;
                    mnb_s = mnb_s/(double)n_cell;
                    
                    Lx1 = Ly1 = ra = rb1 = rb2 = 0;                
                    for (int jj = 0; jj < n_cell; ++jj) {
                        Lx1 += pow(tmpa_s[jj] - mna_s, 2);
                        Ly1 += pow(tmpb_s[jj] - mnb_s, 2);
                        
                        tmpa_s[jj] = tmpa_s[jj] - mna_s;
                        tmpb_s[jj] = tmpb_s[jj] - mnb_s;
                        
                        ra += tmpa_s[jj] * tmpb_s[jj];
                        rb1 += pow(tmpa_s[jj],2);
                        rb2 += pow(tmpb_s[jj],2);
                    }
                    
                    rb1 = sqrt(rb1);
                    rb2 = sqrt(rb2);
                    r = ra/(rb1*rb2);
                    
                    Lx = Lx1/Lx2;
                    Ly = Ly1/Ly2;
                    
                    Ls[k] = sqrt(Lx) * sqrt(Ly) * r;
                    Ds[k] = sqrt(Lx) * (1- r);
                    mn += Ls[k];
                    mn_d += Ds[k];
                }
                mn = mn/perm;
                mn_d = mn_d/perm;
                
                for (int jj = 0; jj < perm; ++jj) {
                    var += pow((Ls[jj]-mn), 2);
                    var_d += pow((Ds[jj]-mn_d), 2);
                }
                var = sqrt(var/perm);
                var_d = sqrt(var_d/perm);
                
                double t = (L-mn)/var;
                double t_d = (D-mn_d)/var_d;
                v->t = t;
                v->L = L;
                v->r = pr;
                v->d = D;
                v->t_d = t_d;
                R_Free(Ls);
                R_Free(Ds);
                R_Free(tmpa_s);
                R_Free(tmpb_s);
                R_Free(tmpa);
                R_Free(tmpb);
            }
            
#pragma omp critical
            {
                dat[i] = vs;
            }
        }
        //
        for (int i = 0; i < n_feature2; ++i) {
            struct val *vs = dat[i];
            for (int j = 0; j < n_feature1; ++j) {
                struct val *v = &vs[j];
                INTEGER(val1)[nn] = v->v1+1;
                INTEGER(val2)[nn] = v->v2+1;
                REAL(val3)[nn] = v->L;
                REAL(val4)[nn] = v->t;
                REAL(val5)[nn] = v->r;
                REAL(val6)[nn] = v->d;
                REAL(val7)[nn] = v->t_d;

                nn++;
            }
            R_Free(vs);
        }
        R_Free(dat);
    }

    for (int pi = 0; pi < perm; ++pi) R_Free(ris[pi]);
    R_Free(ris);
    
    SEXP ta = PROTECT(allocVector(VECSXP, 7));
    SET_VECTOR_ELT(ta, 0, val1);
    SET_VECTOR_ELT(ta, 1, val2);
    SET_VECTOR_ELT(ta, 2, val3);
    SET_VECTOR_ELT(ta, 3, val4);
    SET_VECTOR_ELT(ta, 4, val5);
    SET_VECTOR_ELT(ta, 5, val6);
    SET_VECTOR_ELT(ta, 6, val7);
    UNPROTECT(8);

    return ta;
}
