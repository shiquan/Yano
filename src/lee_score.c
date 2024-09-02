#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>

static cholmod_common c;

extern int *random_idx(const int n);
extern void shuffle(double tmp[], int const idx[], const int n);
extern void smooth_W(double * const a, double *s, const int N, CHM_SP W);

struct val {
    int v1; // feature 
    int v2; // test.feature
    double L;
    double t;
};
SEXP lee_score(SEXP _A, SEXP _B, SEXP _W, SEXP _perm, SEXP _threads, SEXP _seed)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP B = AS_CHM_SP__(_B);
    CHM_SP W = AS_CHM_SP__(_W);

    const int perm = asInteger(_perm);
    int n_thread = asInteger(_threads);

    const int seed = asInteger(_seed);

    srand(seed);
    
    if (A->ncol != B->ncol) return mkString("Unequal column of A and B.");
    if (A->ncol != W->nrow) return mkString("A column and W row do not match.");
    if (W->ncol != W->nrow) return mkString("W is not a square matrix.");

    const int *ap = (int *)A->p;
    const int *ai = (int *)A->i;
    const double *ax = (double *)A->x;

    const int *bp = (int *)B->p;
    const int *bi = (int *)B->i;
    const double *bx = (double *)B->x;

    const int n_cell = A->ncol;
    
    int **ris = NULL;
    ris = R_Calloc(perm, int*);
    for (int pi = 0; pi < perm; ++pi) {
        ris[pi] = random_idx(n_cell);
    }

    R_CheckStack();

    R_CheckUserInterrupt();

    int nn = 0;
    
    const int n_feature1 = A->nrow;
    const int n_feature2 = B->nrow;

    int n_record = n_feature1 * n_feature2;
    
    SEXP val1 = PROTECT(allocVector(INTSXP, n_record)); // feature
    SEXP val2 = PROTECT(allocVector(INTSXP, n_record)); // test.feature
    SEXP val3 = PROTECT(allocVector(REALSXP, n_record)); // L
    SEXP val4 = PROTECT(allocVector(REALSXP, n_record)); // t
    
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
                for (int k = 0; k < n_cell; ++k) {
                    Lx1 += pow(tmpa_s[k]-mna,2);
                    Lx2 += pow(tmpa[k]-mna,2);
                    Ly1 += pow(tmpb_s[k]-mnb,2);
                    Ly2 += pow(tmpb[k]-mnb,2);
                    
                    tmpa_s[k] = tmpa_s[k] - mna_s;
                    tmpb_s[k] = tmpb_s[k] - mnb_s;
                    ra  += tmpa_s[k] * tmpb_s[k];
                    ra  += tmpa_s[k] * tmpb[k];
                    rb1 += pow(tmpa_s[k],2);
                    rb2 += pow(tmpb_s[k],2);
                }
                
                rb1 = sqrt(rb1);
                rb2 = sqrt(rb2);
                
                double Lx = Lx1/Lx2;
                double Ly = Ly1/Ly2;
                double r = ra/(rb1*rb2);
                double L = sqrt(Lx) * sqrt(Ly) *r ;

                double mn = 0, var = 0;
                double *Ls = R_Calloc(perm, double);
                
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
                    mn += Ls[k];
                }
                mn = mn/perm;
                
                for (int jj = 0; jj < perm; ++jj) {
                    var += pow((Ls[jj]-mn), 2);
                }
                var = sqrt(var/perm);
                
                double t = (L-mn)/var;
                v->t = t;
                v->L = L;
                // ts[j] = t;
                // LL[j] = L;
                
                R_Free(Ls);
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
        for (int i = 0; i < n_feature1; ++i) {
            struct val *vs = dat[i];
            for (int j = 0; j < n_feature2; ++j) {
                struct val *v = &vs[j];
                INTEGER(val1)[nn] = v->v1+1;
                INTEGER(val2)[nn] = v->v2+1;
                REAL(val3)[nn] = v->L;
                REAL(val4)[nn] = v->t;
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
                mna_s = mna_s/(double)n_cell;
                mnb_s = mnb_s/(double)n_cell;
                
                double Lx1 = 0, Lx2 = 0, Ly1 = 0, Ly2 = 0, ra = 0, rb1 = 0, rb2 = 0;
                for (int k = 0; k < n_cell; ++k) {
                    Lx1 += pow(tmpa_s[k]-mna,2);
                    Lx2 += pow(tmpa[k]-mna,2);
                    Ly1 += pow(tmpb_s[k]-mnb,2);
                    Ly2 += pow(tmpb[k]-mnb,2);
                    
                    tmpa_s[k] = tmpa_s[k] - mna_s;
                    tmpb_s[k] = tmpb_s[k] - mnb_s;
                    ra  += tmpa_s[k] * tmpb_s[k];
                    ra  += tmpa_s[k] * tmpb[k];
                    rb1 += pow(tmpa_s[k],2);
                    rb2 += pow(tmpb_s[k],2);
                }
                
                rb1 = sqrt(rb1);
                rb2 = sqrt(rb2);
                
                double Lx = Lx1/Lx2;
                double Ly = Ly1/Ly2;
                double r = ra/(rb1*rb2);
                double L = sqrt(Lx) * sqrt(Ly) *r ;

                double mn = 0, var = 0;
                double *Ls = R_Calloc(perm, double);
                
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
                    mn += Ls[k];
                }
                mn = mn/perm;
                
                for (int jj = 0; jj < perm; ++jj) {
                    var += pow((Ls[jj]-mn), 2);
                }
                var = sqrt(var/perm);
                
                double t = (L-mn)/var;
                v->t = t;
                v->L = L;
                
                R_Free(Ls);
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
                nn++;
            }
            R_Free(vs);
        }
        R_Free(dat);
    }

    for (int pi = 0; pi < perm; ++pi) R_Free(ris[pi]);
    R_Free(ris);
    
    SEXP ta = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ta, 0, val1);
    SET_VECTOR_ELT(ta, 1, val2);
    SET_VECTOR_ELT(ta, 2, val3);
    SET_VECTOR_ELT(ta, 3, val4);

    UNPROTECT(5);

    return ta;
}
