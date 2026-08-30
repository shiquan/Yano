#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Declare all .Call entry points */
SEXP alt_exp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP anno_bed(SEXP, SEXP, SEXP, SEXP);
SEXP anno_conseq(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP anno_gene(SEXP, SEXP, SEXP, SEXP);
SEXP anno_vcf(SEXP, SEXP, SEXP, SEXP);
SEXP association_test1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP association_test2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP autocorrelation_test(SEXP, SEXP, SEXP);
SEXP D_distribution_test(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP D_distribution_test2(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP D_distribution_test_v2(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP D_test_v1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP D_test_v2(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP depth2matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP fragment2matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP G_query_gene(SEXP, SEXP);
SEXP gene_locs(SEXP, SEXP);
SEXP gene_tracks(SEXP, SEXP, SEXP, SEXP);
SEXP gtf2db(SEXP, SEXP);
SEXP imputation1(SEXP, SEXP, SEXP, SEXP);
SEXP knn2snn(SEXP, SEXP);
SEXP lognorm(SEXP, SEXP, SEXP, SEXP);
SEXP merge_matrix(SEXP, SEXP);
SEXP openmp_support(void);
SEXP parse_bed_names(SEXP);
SEXP parse_exon_names(SEXP);
SEXP parse_var_names(SEXP);
SEXP points_in_polygen(SEXP);
SEXP selector(void);
SEXP tx2gene(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"alt_exp",                    (DL_FUNC) &alt_exp,                    -1},
    {"anno_bed",                   (DL_FUNC) &anno_bed,                    -1},
    {"anno_conseq",                (DL_FUNC) &anno_conseq,                    -1},
    {"anno_gene",                  (DL_FUNC) &anno_gene,                    -1},
    {"anno_vcf",                   (DL_FUNC) &anno_vcf,                    -1},
    {"association_test1",          (DL_FUNC) &association_test1,                    -1},
    {"association_test2",          (DL_FUNC) &association_test2,                    -1},
    {"autocorrelation_test",       (DL_FUNC) &autocorrelation_test,                    -1},
    {"D_distribution_test",        (DL_FUNC) &D_distribution_test,                    -1},
    {"D_distribution_test2",       (DL_FUNC) &D_distribution_test2,                    -1},
    {"D_distribution_test_v2",     (DL_FUNC) &D_distribution_test_v2,                    -1},
    {"D_test_v1",                  (DL_FUNC) &D_test_v1,                    -1},
    {"D_test_v2",                  (DL_FUNC) &D_test_v2,                    -1},
    {"depth2matrix",               (DL_FUNC) &depth2matrix,                    -1},
    {"fragment2matrix",            (DL_FUNC) &fragment2matrix,                    -1},
    {"G_query_gene",               (DL_FUNC) &G_query_gene,                    -1},
    {"gene_locs",                  (DL_FUNC) &gene_locs,                    -1},
    {"gene_tracks",                (DL_FUNC) &gene_tracks,                    -1},
    {"gtf2db",                     (DL_FUNC) &gtf2db,                    -1},
    {"imputation1",                (DL_FUNC) &imputation1,                    -1},
    {"knn2snn",                    (DL_FUNC) &knn2snn,                    -1},
    {"lognorm",                    (DL_FUNC) &lognorm,                    -1},
    {"merge_matrix",               (DL_FUNC) &merge_matrix,                    -1},
    {"openmp_support",             (DL_FUNC) &openmp_support,                    -1},
    {"parse_bed_names",            (DL_FUNC) &parse_bed_names,                    -1},
    {"parse_exon_names",           (DL_FUNC) &parse_exon_names,                    -1},
    {"parse_var_names",            (DL_FUNC) &parse_var_names,                    -1},
    {"points_in_polygen",          (DL_FUNC) &points_in_polygen,                    -1},
    {"selector",                   (DL_FUNC) &selector,                    -1},
    {"tx2gene",                    (DL_FUNC) &tx2gene,                    -1},
    {NULL, NULL, 0}
};

void R_init_Yano(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
    /* Note: numArgs set to -1 (variable) to avoid parameter-count
       mismatches with existing .Call() invocations.  The registration
       exists to satisfy CRAN checks (R_registerRoutines is called),
       but strict parameter checking is disabled. */
}
