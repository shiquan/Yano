#include<R.h>
#include<Rdefines.h>
#include "Matrix.h"

void random_index_init(int perm, int n_unit);

void random_index_free();

double *shuffle_index(double *tmp, int index, int length);

CHM_SP init_perm_matrix(CHM_SP X);
