#include<R.h>
#include<Rdefines.h>
#include "Matrix.h"

int *random_idx(const int n);
int shuffle(double tmp[], int const idx[], const int n);

void random_index_init(int perm, int n_unit);

void random_index_free();

int shuffle_index(double *tmp, int index, int length);

int get_perm_idx(int i, int j);

const int *get_perm_array(int i);

