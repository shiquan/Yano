#include<R.h>
#include<Rdefines.h>

void random_index_init(int perm, int n_unit);

void random_index_free();

double *shuffle_index(double *tmp, int index, int length);
