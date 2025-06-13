#include "perm.h"

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

static int **ris = NULL;
static int permutation = 0;

void random_index_init(int perm, int n_unit) {
    permutation = perm;
    ris = R_Calloc(perm,int*);
    for (int pi = 0; pi < perm; ++pi) {
        ris[pi] = random_idx(n_unit);
    }    
}

void random_index_free()
{
    int i;
    for (i = 0; i < permutation; ++i) R_Free(ris[i]);
    R_Free(ris);
    permutation = 0;
    ris = NULL;
}

double *shuffle_index(double *tmp, int index, int length)
{
    shuffle(tmp, ris[index], length);
    return tmp;
}

CHM_SP init_perm_matrix(CHM_SP X)
{
    static cholmod_common c;
    
    int *xi = (int*)X->i;
    int *xp = (int*)X->p;
    int i, p;
    for (i = 1; i < permutation+1; ++i) { // start from 1, because 0 is orginal X
        int *idx = ris[i];
        for (p = xp[i]; p < xp[i+1]; ++p) {
            xi[p] = idx[xi[p]];
        }
    }

    M_cholmod_sort(X, &c);
}
