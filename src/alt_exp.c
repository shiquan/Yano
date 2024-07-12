#ifdef _OPENMP
#include <omp.h>
#endif

#include<R.h>
#include<Rdefines.h>
#include <Matrix.h>
// # nocov end
static cholmod_common c;

extern int *random_idx(const int n);
extern void shuffle(double tmp[], int const idx[], const int n);

SEXP alt_exp(SEXP _A, SEXP _B, SEXP idx1, SEXP idx2, SEXP _mode, SEXP _perm, SEXP _threads, SEXP _seed)
{
    CHM_SP A = AS_CHM_SP__(_A);
    CHM_SP B = AS_CHM_SP__(_B);

    const int perm = asInteger(_perm);
    const int n_thread = asInteger(_threads);

    const int seed = asInteger(_seed);

    srand(seed);
    
    int mode = asInteger(_mode);
    if (A->stype) return mkString("A cannot be symmetric");
    if (B->stype) return mkString("B cannot be symmetric");
    
    if (A->ncol != B->ncol) return mkString("A column and B column do not match.");
    if (A->nrow != B->nrow) return mkString("A row and B row do not match.");

    if (A->ncol < 2) return mkString("Too few cells."); // to do

    A = M_cholmod_transpose(A, (int)A->xtype, &c);
    B = M_cholmod_transpose(B, (int)B->xtype, &c);

    R_CheckStack();

    const int n_cell = A->nrow;
    const int N_feature = A->ncol;

    const int len_g1 = length(idx1);
    const int len_g2 = length(idx2);

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
        
#pragma omp parallel for num_threads(n_thread)
    for (int i = 0; i < N_feature; ++i) {
        double *tmpa = R_Calloc(n_cell, double);
        double *tmpb = R_Calloc(n_cell, double);
        memset(tmpa, 0, sizeof(double)*n_cell);
        memset(tmpb, 0, sizeof(double)*n_cell);

        for (int j = ap[i]; j < ap[i+1]; ++j) {
            if (ISNAN(ax[j])) continue;
            int cid = ai[j];
            tmpa[cid] = ax[j];
        }
        
        for (int j = bp[i]; j < bp[i+1]; ++j) {
            if (ISNAN(bx[j])) continue;
            int cid = bi[j];
            tmpb[cid] = bx[j];
        }

        double g11 = 0;
        double g21 = 0;
        double g12 = 0;
        double g22 = 0;

        for (int j = 0; j < len_g1; ++j) {
            int ii = INTEGER(idx1)[j]  -1;
            g11 += tmpa[ii];
            g12 += tmpb[ii];
        }
        for (int j = 0; j < len_g2; ++j) {
            int ii = INTEGER(idx2)[j]  -1;
            g21 += tmpa[ii];
            g22 += tmpb[ii];
        }

        if (mode == 2) {
            g12 = g12 - g11;
            g22 = g22 - g21;
        } else if (mode == 3) {
            g12 = g12 + g11;
            g22 = g22 + g21;
        }

        //if (g12 == 0) g12 = 0.001;
        //if (g22 == 0) g22 = 0.001;

        double delta = (g11+1)/(g12+1) - (g21+1)/(g22+1);
        
        double mean = 0;
        double var = 0;
        
        double *es = R_Calloc(perm, double);
        for (int k = 0; k < perm; ++k) {
            shuffle(tmpa, ris[k], n_cell);
            shuffle(tmpb, ris[k], n_cell);

            for (int j = 0; j < len_g1; ++j) {
                int ii = INTEGER(idx1)[j]  -1;
                g11 += tmpa[ii];
                g12 += tmpb[ii];
            }
            for (int j = 0; j < len_g2; ++j) {
                int ii = INTEGER(idx2)[j]  -1;
                g21 += tmpa[ii];
                g22 += tmpb[ii];
            }
            
            if (mode == 2) {
                g12 = g12 - g11;
                g22 = g22 - g21;
            } else if (mode == 3) {
                g12 = g12 + g11;
                g22 = g22 + g21;
            }
            
            es[k] = (g11+1)/(g12+1) - (g21+1)/(g22+1);
            mean = mean + es[k];
        }
        R_Free(tmpa);
        R_Free(tmpb);
        
        mean = mean/perm;
        
        for (int k = 0; k < perm; ++k) {
            var += pow((es[k] -mean),2);
        }
        var = sqrt(var/perm);

        double t = (delta - mean)/var;

        R_Free(es);
        
#pragma omp critical
        {
            REAL(Dval)[i] = delta;
            REAL(Tval)[i] = t;
            REAL(Mval)[i] = mean;
            REAL(Vval)[i] = var;
        }
    }

    for (int pi = 0; pi < perm; ++pi) R_Free(ris[pi]);
    R_Free(ris);

    M_cholmod_free_sparse(&A, &c);
    M_cholmod_free_sparse(&B, &c);
    
    SEXP ta = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(ta, 0, Dval);
    SET_VECTOR_ELT(ta, 1, Tval);
    SET_VECTOR_ELT(ta, 2, Mval);
    SET_VECTOR_ELT(ta, 3, Vval);

    UNPROTECT(5);
    return ta;

}
