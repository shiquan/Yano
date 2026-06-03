#include "sparse_ops.h"
#include <stdlib.h>
#include <string.h>

/* ---- internal helpers ---- */

static int ensure_cap(struct sparse_vec *v, int need)
{
    if (v->n_alloc >= need) return 0;
    int newcap = v->n_alloc ? v->n_alloc * 2 : 64;
    if (newcap < need) newcap = need;
    int    *new_idx = realloc(v->idx,  newcap * sizeof(int));
    double *new_val = realloc(v->val, newcap * sizeof(double));
    if (!new_idx || !new_val) {
        free(new_idx); free(new_val);
        return 1;
    }
    v->idx = new_idx;
    v->val = new_val;
    v->n_alloc = newcap;
    return 0;
}

/* ---- public API ---- */

int sparse_vec_from_csc(struct sparse_vec *v,
                        const int *ap, const int *ai, const double *ax,
                        int col)
{
    int nz = ap[col + 1] - ap[col];
    if (ensure_cap(v, nz)) return 1;
    v->nz = 0;
    for (int p = ap[col]; p < ap[col + 1]; p++) {
        if (ISNAN(ax[p])) continue;
        v->idx[v->nz] = ai[p];
        v->val[v->nz] = ax[p];
        v->nz++;
    }
    return 0;
}

int sparse_vec_copy(struct sparse_vec *dst, const struct sparse_vec *src)
{
    if (ensure_cap(dst, src->nz)) return 1;
    memcpy(dst->idx, src->idx, src->nz * sizeof(int));
    memcpy(dst->val, src->val, src->nz * sizeof(double));
    dst->nz = src->nz;
    return 0;
}

void sparse_vec_shuffle(struct sparse_vec *v, const int *perm)
{
    for (int k = 0; k < v->nz; k++)
        v->idx[k] = perm[v->idx[k]];
}

void sparse_smooth(const struct sparse_vec *v,
                   const int *wp, const int *wi, const double *wx,
                   double *result, int n_cell, double filter)
{
    memset(result, 0, n_cell * sizeof(double));

    /* Only iterate non-zero source cells.  For scRNA-seq where most
       features are expressed in &lt;10% of cells this is the dominant
       performance advantage over a dense loop. */
    for (int k = 0; k < v->nz; k++) {
        int    cid = v->idx[k];
        double val = v->val[k];
        for (int j = wp[cid]; j < wp[cid + 1]; j++)
            result[wi[j]] += wx[j] * val;
    }

    /* Apply filter: small values → 0 */
    if (filter > 0.0) {
        for (int i = 0; i < n_cell; i++)
            if (result[i] < filter)
                result[i] = 0.0;
    }
}

void sparse_vec_destroy(struct sparse_vec *v)
{
    free(v->idx);  v->idx = NULL;
    free(v->val);  v->val = NULL;
    v->nz = v->n_alloc = 0;
}
