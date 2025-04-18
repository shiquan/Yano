% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vis.R
\name{PSIPlot}
\alias{PSIPlot}
\title{PSIPlot}
\usage{
PSIPlot(
  object = NULL,
  exon.assay = NULL,
  exclude.assay = "exclude",
  features = NULL,
  dims = c(1, 2),
  cells = NULL,
  cols = c("lightgrey", "red"),
  pt.size = NULL,
  alpha = 1,
  order = TRUE,
  reduction = NULL,
  shape.by = NULL,
  ncol = NULL,
  split.by = NULL,
  by.col = TRUE,
  coord.fixed = FALSE,
  combine = TRUE,
  raster = NULL,
  raster.dpi = c(512, 512)
)
}
\arguments{
\item{object}{Seurat object.}

\item{exon.assay}{Exon assay name.}

\item{exclude.assay}{Excluded exon assay name.}

\item{features}{Features to plot.}

\item{dims}{Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions}

\item{cells}{Vector of cells to plot (default is all cells)}

\item{cols}{Vector of colors, each color corresponds to an identity class. This may also be a single character or numeric value corresponding to a palette as specified by RColorBrewer::brewer.pal.info.}

\item{pt.size}{Adjust point size for plotting}

\item{alpha}{Alpha value for points}

\item{order}{Boolean determing whether to plot cells in order of PSI score.}

\item{reduction}{Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca}

\item{shape.by}{If NULL, all points are circles (default). You can specify any cell attribute (that can be pulled with FetchData) allowing for both different colors and different shapes on cells.  Only applicable if \code{raster = FALSE}.}

\item{ncol}{Number of columns to combine multiple features plots to}

\item{split.by}{A factor in object metadata to split the plot by, pass 'ident' to split by cell identity}

\item{by.col}{If splitting by a factor, plot the splits per column with the features as rows}

\item{coord.fixed}{Plot cartesian coordinates with fixed aspect ratio}

\item{combine}{Combine plots into a single patchwork ggplot object. If \code{FALSE}, return a list of ggplot objects.}

\item{raster.dpi}{Pixel resolution for rasterized plots, passed to geom_scattermore(). Default is c(512, 512).}

\item{min.cutoff, max.cutoff}{Vector of minimum and maximum cutoff values for each feature}

\item{mode}{Test mode. For mode 1, X (test feature) vs Y (binding feature). For mode 2, X vs (Y-X). For mode 3, X vs (Y+X). Please note, when set to mode 2 or 3, will use raw counts to update expression value of binding features. Then normalise the counts before testing. For mode 1, will use Layer 'data'. Default is mode 1.}

\item{group.by}{Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class}
}
\value{
A patchwork ggplot object of \code{combine = TRUE}; otherwise, a list of ggplot objects
}
\description{
Plot PSI score on reduction map.
}
