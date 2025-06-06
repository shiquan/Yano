% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{buildKNN}
\alias{buildKNN}
\title{buildKNN}
\usage{
buildKNN(
  data,
  query = data,
  k.param = 20,
  metric = c("euclidean", "cosine", "manhattan", "hamming"),
  n.trees = 50
)
}
\arguments{
\item{k.param}{K-nearest neighbors, for calculating weight matrix with emb. Default is 20.}

\item{n.trees}{Number of trees during index build time. More trees gives higher accuracy}
}
\description{
Compute KNN with RcppAnnoy.
}
