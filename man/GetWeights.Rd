% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{GetWeights}
\alias{GetWeights}
\title{GetWeights}
\usage{
GetWeights(
  snn = NULL,
  pos = NULL,
  order.cells = NULL,
  emb = NULL,
  k.nn = 20,
  prune.distance = 20,
  prune.SNN = 1/15,
  diag.value = 0,
  cells = NULL
)
}
\arguments{
\item{snn}{Shared nearest neighbour graph, usually can found at object[["RNA_snn"]]. This graph can be calculate by Seurat::FindNeighbors().}

\item{pos}{Tissue coordinates matrix.}

\item{order.cells}{Predefined cell ranks, used for cell lineage analysis.}

\item{emb}{Cell dimesional space (PCA/ICA/harmony).}

\item{k.nn}{K-nearest neighbors, for calculating weight matrix with emb.}

\item{prune.distance}{Sets the cutoff for cell distance on lineage trajectory (ranked cells) or spatial cooridates (bin distance) when computing the neighborhood overlap for the weight matrix construction. Any edges with values greater than this will be set to 0 and removed from the weight matrix graph. Default is 20, means only calculate weight edges for nearby 20 cells for each cell.}

\item{prune.SNN}{Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 --- prune everything). Default is 1/15.}

\item{diag.value}{Diagnoal value in the weight matrix.}

\item{cells}{Cell list. Default use all cells.}
}
\value{
A sparse weight matrix.
}
\description{
Calcualte cell-cell weight matrix by one of shared nearest neighbour matrix, spatial locations, cell embedding and linear trajectory.
}
