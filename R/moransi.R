#' Calculate Moran's index for features in parallel.
#' @param object
#' @param assay
#' @param slot
#' @param spatial
#' @param W
#' @param scale.weight
#' @param reduction
#' @param dims
#' @param prune.distance Sets cutoff for distance matrix. Any edges with distance
#' greater than this value will be set to 0 and remove from the distance matrix.
#' @param order.cells
#' @param cells
#' @param features
#' @param threads
#' @param ... parameters pass to
#' @importFrom methods as
#' @export
RunAutoCorr <- function(object = NULL,
                        assay = NULL,                       
                        slot = "data",
                        spatial = FALSE,
                        snn.name = NULL,
                        order.cells = NULL,
                        prune.distance = 20,
                        prune.snn = 1/10,
                        cells = NULL,
                        features = NULL,
                        weight.matrix.name = "WeightMatrix",
                        prefix="moransi",
                        #perm = 100,
                        threads=0,
                        ...)
{
  check.par <- 0
  if (!is.null(snn.name)) check.par <- check.par + 1
  if (!is.null(order.cells)) check.par <- check.par + 1
  if (isTRUE(spatial)) check.par <- check.par + 1

  if (check.par > 1) {
    stop("Only support to set one of spatial, snn.name, and order.cells.")
  }
  
  if (check.par == 0) snn.name <- "RNA_snn"

  tt <- Sys.time()
  
  assay <- assay %||% DefaultAssay(object = object)
  message(paste0("Working on assay : ", assay))

  cells <- order.cells %||% cells
  cells <- cells %||% colnames(object)
  cells <- intersect(colnames(object), cells)
  
  features <- features %||% rownames(object)
  features <- intersect(rownames(object),features)

  threads <- getCores(threads)

  if (!is.null(snn.name)) {
    W <- GetWeightsFromSNN(object=object, name = snn.name, prune.snn = prune.snn)
  }

  if (!is.null(order.cells)) {
    W <- GetWeights(order.cells = order.cells, prune.distance = prune.distance)
  }

  if (isTRUE(spatial)) {
    W <- GetWeightsFromSpatial(object = object, prune.distance = prune.distance, ...)
  }

  W <- as(W, "Graph")
  object[[weight.matrix.name]] <- W
  
  x0 <- GetAssayData(object, assay = assay, slot = slot)[,cells]
  features <- intersect(rownames(x0), features)
  x0 <- x0[features,]
  x0 <- as(x0, "dgCMatrix")
  ## if (perm < 10) {
  ##   perm <- 0
  ## }
  
  message(paste0("Run autocorrelation test for ", length(features), " features."))
  #moransi.vals <- .Call("moransi_perm_test", x0, W, TRUE, threads, FALSE, perm)
  moransi.vals <- .Call("autocorrelation_test", x0, W, TRUE, threads)

  if (length(moransi.vals) == 1) stop(moransi.vals[[1]])

  Ivals <- moransi.vals[[1]]
  names(Ivals) <- features

  ## if (perm > 1) {
  ##   Tvals <- moransi.vals[[2]]
  ##   names(Tvals) <- features
  ##   pvals <- pt(Tvals, df = perm - 1, lower.tail = FALSE)
  ##   names(pvals) <- features
  ##   prefix.p <- paste0(prefix, ".pval")
  ##   object[[assay]]@meta.features[[prefix.p]] <- pvals[rownames(object)]
  ## }

  Zvals <- moransi.vals[[2]]
  names(Zvals) <- features
  pvals <- pnorm(Zvals, lower.tail = FALSE)
  names(pvals) <- features
  prefix.p <- paste0(prefix, ".pval")
  object[[assay]]@meta.features[[prefix.p]] <- pvals[rownames(object)]
  object[[assay]]@meta.features[[prefix]] <- Ivals[rownames(object)]
  
  rm(moransi.vals)
  gc()

  tt <- Sys.time()-tt
  message(paste0("Runtime : ",format(tt)));
  object
}

#' Set autocorrection by rank and/or score
#' @export
SetAutoCorrFeatures <- function(object = NULL,
                                moransi.min = 0,
                                assay = DefaultAssay(object),
                                p.cutoff = 0.01,
                                prefix = "moransi")
{
  tab <- object[[assay]]@meta.features

  cn <- colnames(tab)
  prefix.p <- paste0(prefix,".pval")
  if (prefix %ni% cn | prefix.p %ni% cn) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }

  idx <- which(tab[[prefix]] > moransi.min & tab[[prefix.p]] <= p.cutoff)

  message(paste0(length(idx), " autocorrelated features."))
  tab[["autocorr.variable"]] <- FALSE
  tab[idx,][["autocorr.variable"]] <- TRUE
  object[[assay]]@meta.features <- tab
  
  object
}
#' @export
AutoCorrFeatures <- function(object = NULL, assay = NULL)
{
  assay <- assay %||% DefaultAssay(object)
  tab <- object[[assay]]@meta.features

  if ('autocorr.variable' %ni% colnames(tab)) {
    stop("No autocorrelation flag found, run SetAutoCorrFeatures() first.")
  }
  rownames(object)[which(tab[['autocorr.variable']])]
}
