% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dissimilarity.R
\name{RunBlockCorr}
\alias{RunBlockCorr}
\title{Run spatial dissimilarity test for features and their binding features in parallel.}
\usage{
RunBlockCorr(
  object = NULL,
  bind.name = "gene_name",
  features = NULL,
  assay = NULL,
  min.cells = 10,
  bind.assay = NULL,
  bind.features = NULL,
  min.cells.bind = 10,
  prefix = NULL,
  subset = NULL,
  min.features.per.block = 1,
  scale.factor = 10000,
  weight.matrix.name = "WeightMatrix",
  mode = c(1, 2, 3),
  method = c("D1", "D2"),
  library.size = NULL,
  scale = FALSE,
  perm = 100,
  seed = 999,
  threads = 0,
  block.name = NULL,
  block.assay = NULL,
  block.features = NULL,
  verbose = TRUE,
  debug = FALSE
)
}
\arguments{
\item{object}{Seurat object}

\item{bind.name}{Title name for binding features in the meta table. Consider most users start Yano to perform alternative splicing analysis, the default bind.name set to "gene_name".}

\item{features}{Vector of features to calculate. Default is AutoCorrFeatures(object).}

\item{assay}{Work assay.}

\item{min.cells}{Features detected in few than minimal number of cells will be skipped. Default is 10.}

\item{bind.assay}{Name of binding assay.}

\item{bind.features}{List of bind features. Default use all.}

\item{min.cells.bind}{Binding features detected in few than minimal number of cells will be skipped. Default is 10.}

\item{prefix}{Prefix name for output scores and values. Default is same with bind.name.}

\item{subset}{Rules for subset meta table before select features to perform test.}

\item{scale.factor}{Scale factor to normalise counts. Default is 1e4. For mode 1, this function will use data from Layer 'data'. For mode 2 and 3, will use data from Layer 'counts'. The counts will further be normalised with sample size and scale factor for spatial dissimilarity test.}

\item{weight.matrix.name}{Weight matrix name, this matrix (graph) generated by \code{\link{RunAutoCorr}}.}

\item{mode}{Test mode. For mode 1, X (test feature) vs Y (binding feature). For mode 2, X vs (Y-X). For mode 3, X vs (Y+X). Please note, when set to mode 2 or 3, will use raw counts to update expression value of binding features. Then normalise the counts before testing. For mode 1, will use Layer 'data'. Default is mode 1.}

\item{method}{Method to use. \eqn{D1 = sqrt(L_x)(1-r_xy)}, \eqn{D2 = sqrt(L_x)sqrt(L_y)(1-r_xy)}. In default use D1 method, see the manual for details.}

\item{library.size}{Library size for each cell, used for normalise counts when mode is 2 or 3. If not set, use colSum(counts) instead.}

\item{scale}{Scale expression of the test feature and binding feature when perfroming test.}

\item{perm}{Permutations for evaluating mean and sd of D scores. Default is 100.}

\item{seed}{Seed for generate random number. Default is 999.}

\item{threads}{Threads. If set to 0 (default), will auto check the CPU cores and set threads = number of CPU cores -1.}

\item{versbose}{Print log message. Default is TRUE.}
}
\description{
Run spatial dissimilarity test for features and their binding features in parallel.
}
