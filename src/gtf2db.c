#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "dict.h"
#include "gtf.h"
#include "bed.h"

static void
_finalizer(SEXP ext)
{
    struct gtf_spec *G = (struct gtf_spec*)R_ExternalPtrAddr(ext);
    anno_bed_cleanup();
    gtf_destroy(G);
}

SEXP gtf2db(SEXP filename, SEXP utr)
{
    if (!Rf_isString(filename)) error_return("Input is not String");
    const char *file = translateChar(STRING_ELT(filename, 0));

    int use_utr = asInteger(utr);
    fprintf(stderr, "%s\n", file);
    struct gtf_spec *G = gtf_read(file, use_utr?1:2);

    SEXP ext = PROTECT(R_MakeExternalPtr(G, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ext, _finalizer, TRUE);
    UNPROTECT(1);
    return ext;
}

