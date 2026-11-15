#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Declare all .Call entry points.  Argument counts must match the
   definitions in the corresponding source files. */
SEXP alt_exp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP anno_bed(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_anno_conseq(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_anno_gene(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_anno_vcf(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP association_test1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP association_test2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP autocorrelation_test(SEXP, SEXP, SEXP, SEXP);
SEXP D_distribution_test(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP D_distribution_test2(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP D_distribution_test_v2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP D_test_v1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP D_test_v2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP depth2matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP fragment2matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP G_query_gene(SEXP, SEXP);
SEXP gene_locs(SEXP, SEXP);
SEXP gene_tracks(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP C_gtf2db(SEXP, SEXP);
SEXP imputation1(SEXP, SEXP, SEXP, SEXP);
SEXP knn2snn(SEXP, SEXP);
SEXP C_lognorm(SEXP, SEXP, SEXP);
SEXP merge_matrix(SEXP);
SEXP openmp_support(void);
SEXP parse_bed_names(SEXP);
SEXP parse_exon_names(SEXP);
SEXP parse_var_names(SEXP);
SEXP points_in_polygen(SEXP, SEXP, SEXP, SEXP);
SEXP selector(void);
SEXP tx2gene(SEXP, SEXP);

/* Registered routine names are the objects that useDynLib assigns into
   the package namespace.  The five annotation/normalisation routines are
   prefixed with C_ so they cannot collide with the same-named R wrappers
   (anno_vcf, anno_gene, anno_conseq, gtf2db, lognorm). */
static const R_CallMethodDef CallEntries[] = {
    {"alt_exp",                    (DL_FUNC) &alt_exp,                    7},
    {"anno_bed",                   (DL_FUNC) &anno_bed,                   7},
    {"C_anno_conseq",              (DL_FUNC) &C_anno_conseq,              8},
    {"C_anno_gene",                (DL_FUNC) &C_anno_gene,                7},
    {"C_anno_vcf",                 (DL_FUNC) &C_anno_vcf,                 9},
    {"association_test1",          (DL_FUNC) &association_test1,          6},
    {"association_test2",          (DL_FUNC) &association_test2,          6},
    {"autocorrelation_test",       (DL_FUNC) &autocorrelation_test,       4},
    {"D_distribution_test",        (DL_FUNC) &D_distribution_test,        5},
    {"D_distribution_test2",       (DL_FUNC) &D_distribution_test2,       5},
    {"D_distribution_test_v2",     (DL_FUNC) &D_distribution_test_v2,     6},
    {"D_test_v1",                  (DL_FUNC) &D_test_v1,                 14},
    {"D_test_v2",                  (DL_FUNC) &D_test_v2,                  9},
    {"depth2matrix",               (DL_FUNC) &depth2matrix,              14},
    {"fragment2matrix",            (DL_FUNC) &fragment2matrix,            8},
    {"G_query_gene",               (DL_FUNC) &G_query_gene,               2},
    {"gene_locs",                  (DL_FUNC) &gene_locs,                  2},
    {"gene_tracks",                (DL_FUNC) &gene_tracks,                5},
    {"C_gtf2db",                   (DL_FUNC) &C_gtf2db,                   2},
    {"imputation1",                (DL_FUNC) &imputation1,                4},
    {"knn2snn",                    (DL_FUNC) &knn2snn,                    2},
    {"C_lognorm",                  (DL_FUNC) &C_lognorm,                  3},
    {"merge_matrix",               (DL_FUNC) &merge_matrix,               1},
    {"openmp_support",             (DL_FUNC) &openmp_support,             0},
    {"parse_bed_names",            (DL_FUNC) &parse_bed_names,            1},
    {"parse_exon_names",           (DL_FUNC) &parse_exon_names,           1},
    {"parse_var_names",            (DL_FUNC) &parse_var_names,            1},
    {"points_in_polygen",          (DL_FUNC) &points_in_polygen,          4},
    {"selector",                   (DL_FUNC) &selector,                   0},
    {"tx2gene",                    (DL_FUNC) &tx2gene,                    2},
    {NULL, NULL, 0}
};

void R_init_Yano(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
