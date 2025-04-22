#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "utils.h"
#include "number.h"
#include "htslib/kstring.h"
#include "htslib/bgzf.h"
// #include "Matrix/Matrix.h"
#include <Matrix.h>
// # nocov end
static cholmod_common c;

SEXP readmex(SEXP _file)
{
    const char *file;
    BGZF *fp = bgzf_open(file, "r");
    kstring_t str = {0,0,0};
    int n = 0;
    CHM_SP ans = NULL;
    int *ap;
    int *ai;
    double *ax;

    int lst = -1;
    int k = 0;
    int j = 0;
    for (;;) {
        int rst = bgzf_getline(fp, '\n', &str);
        if (rst < 0) break;

        if (str.s[0] == '%') continue;
        int p;
        int *s = ksplit(&str, 0, &p);
        if (p != 3) {
            Rprintf("malformed format.");
            break;
        }

        int n_feature = str2int(str.s+s[0]);
        int n_barcode = str2int(str.s+s[1]);
        int n_record = str2int(str.s+s[2]);

        if (n == 0) {
            ans = M_cholmod_allocate_sparse(n_feature, n_barcode, n_record, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
            if (ans == NULL) return(mkString("Failed to create sparse matrix"));
            
            ap = (int *)ans->p;
            ai = (int *)ans->i;
            ax = (double *)ans->x;
        } else {
            if (lst == -1) {
                lst = n_feature;
                ap[j++] = k;
            }
            if (lst != n_feature) {
                ap[j++] = k;
            }
            ai[k] = n_barcode;
            ax[k] = n_record;
            
            k++;
        }
        n++;
    }

    UNPROTECT(1);
    
    return M_chm_sparse_to_SEXP(ans, 1, -1, 0, "N", R_NilValue);
}
