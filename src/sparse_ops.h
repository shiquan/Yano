#ifndef SPARSE_OPS_H
#define SPARSE_OPS_H

#include <R.h>
#include <stddef.h>

/* A sparse vector: non-zero entries of one feature across cells.
   Thread-safe — no global state, all allocations are caller-owned. */
struct sparse_vec {
    int    nz;      /* number of non-zero entries            */
    int   *idx;     /* cell index for each entry [0..ncell)  */
    double *val;    /* expression value for each entry       */
    int    n_alloc; /* allocated capacity (>= nz)            */
};

/* Extract column `col` from CSC sparse matrix `A` into a sparse_vec.
   Returns 0 on success, 1 on allocation failure. */
int sparse_vec_from_csc(struct sparse_vec *v,
                        const int *ap, const int *ai, const double *ax,
                        int col);

/* Deep copy src → dst. dst must be initialised (or zeroed).
   Returns 0 on success, 1 on allocation failure. */
int sparse_vec_copy(struct sparse_vec *dst, const struct sparse_vec *src);

/* Shuffle cell indices in-place using a pre-computed permutation.
   perm[idx] gives the new cell index for old cell index idx. */
void sparse_vec_shuffle(struct sparse_vec *v, const int *perm);

/* Smooth a sparse vector through weight matrix W:
     result[i] = sum_{j in N(i)} W[i,j] * x[j]
   where x is the sparse vector (zero for cells not in v).
   result must be pre-allocated with n_cell elements.
   Elements of result below `filter` are set to zero. */
void sparse_smooth(const struct sparse_vec *v,
                   const int *wp, const int *wi, const double *wx,
                   double *result, int n_cell, double filter);

/* Release internal buffers. The struct itself is not freed. */
void sparse_vec_destroy(struct sparse_vec *v);

#endif /* SPARSE_OPS_H */
