#' Calculate Moran's index for features in parallel.
#' @param object
#' @param assay
#' @param layer
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
                        layer = "data",
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

  layer <- Layers(object = object, search = layer)
  if (is.null(layer)) {
    abort("No layer found. Please run NormalizeData or RunTFIDF and retry..")
  }
  
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
    W <- GetWeightsFromSNN(object=object, name = snn.name, prune.snn = prune.snn, cells = cells)
  }

  if (!is.null(order.cells)) {
    W <- GetWeights(order.cells = order.cells, prune.distance = prune.distance)
  }

  if (isTRUE(spatial)) {
    W <- GetWeightsFromSpatial(object = object, prune.distance = prune.distance, ...)
  }

  ncell <- ncol(object)
  cells1 <- colnames(W)
  cells2 <- setdiff(colnames(object), cells1)
  cells1 <- c(cells1, cells2)
  W <- as(W, "TsparseMatrix")
  W <- Matrix::sparseMatrix(i = W@i+1, j = W@j+1, x = W@x, dims = c(ncell,ncell))
  colnames(W) <- cells1
  rownames(W) <- cells1
  W <- W[colnames(object), colnames(object)]
  W <- as(W, "Graph")
  object[[weight.matrix.name]] <- W
  
  x0 <- GetAssayData(object, assay = assay, layer = layer)[,cells]
  features <- intersect(rownames(x0), features)
  x0 <- x0[features,]
  x0 <- as(x0, "CsparseMatrix")

  W <- W[cells, cells]
  message(paste0("Run autocorrelation test for ", length(features), " features."))
  moransi.vals <- .Call("autocorrelation_test", x0, W, TRUE, threads)

  if (length(moransi.vals) == 1) stop(moransi.vals[[1]])

  Ivals <- moransi.vals[[1]]
  names(Ivals) <- features

  Zvals <- moransi.vals[[2]]
  names(Zvals) <- features
  pvals <- pnorm(Zvals, lower.tail = FALSE)
  names(pvals) <- features
  prefix.p <- paste0(prefix, ".pval")

  object0 <- object[[assay]]
  object0[[prefix.p]] <- pvals[rownames(object)]
  object0[[prefix]] <- Ivals[rownames(object)]

  object[[assay]] <- object0
  
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
  object0 <- object[[assay]]

  cn <- colnames(object0[[]])
  prefix.p <- paste0(prefix,".pval")
  
  if (prefix %ni% cn | prefix.p %ni% cn) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }

  idx <- which(object0[[prefix]] > moransi.min & object0[[prefix.p]] <= p.cutoff)

  message(paste0(length(idx), " autocorrelated features."))
  object0[["autocorr.variable"]] <- FALSE

  all <- object0[["autocorr.variable"]][[1]]
  all[idx] <- TRUE
  object0[["autocorr.variable"]] <- all
  object[[assay]] <- object0
  
  object
}
#' @export
AutoCorrFeatures <- function(object = NULL, assay = NULL)
{
  assay <- assay %||% DefaultAssay(object)

  object0 <- object[[assay]]

  if ('autocorr.variable' %ni% colnames(object0[[]])) {
    stop("No autocorrelation flag found, run SetAutoCorrFeatures() first.")
  }
  rownames(object)[which(object0[['autocorr.variable']][[1]])]
}
