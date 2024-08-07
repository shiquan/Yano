% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocess.R
\name{ParseVAR}
\alias{ParseVAR}
\title{ParseVAR}
\usage{
ParseVAR(object = NULL, assay = NULL)
}
\arguments{
\item{object}{Seurat object}

\item{assay}{EAT assay. Default is current actived assay.}
}
\description{
Parse chromosome name, start, ref allele, alternative allele and strand from EAT name. EAT generated by `PISA anno` is formated like chr:pos(ref)>(alt)/[+-].
}
\examples{
data("neuron_small")
DefaultAssay(neuron_small) <- "EAT"
# Check the meta table of features before parsing
head(neuron_small[['EAT']][[]])
neuron_small <- ParseVAR(neuron_small)

# Now see the meta table after parsing
head(neuron_small[['EAT']][[]])
}
