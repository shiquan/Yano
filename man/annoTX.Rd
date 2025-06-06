% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/txanno.R
\name{annoTX}
\alias{annoTX}
\title{annoTX}
\usage{
annoTX(object = NULL, assay = NULL, gtf = NULL, gene.name = "gene_name")
}
\arguments{
\item{object}{Seurat object.}

\item{assay}{Work assay.}

\item{gtf}{GTF object, load by gtf2db.}

\item{gene.name}{Tag name for gene name. Will be created after annotation. Default is "gene_name".}
}
\value{
Annotated Seurat object with annotations.
}
\description{
Annotate gene name and genomic locations for transcript id.
}
