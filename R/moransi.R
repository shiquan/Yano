#' @name RunAutoCorr
#' @title Calculate spatial autocorrelation (Moran's index) for features in parallel.
#' @param object Seurat object
#' @param assay Working assay
#' @param layer Input data. For expression, log scaled data is highly recommended
#' @param spatial Use spatial coordinate instead of SNN space and linear trajectory to calculate the cell-cell weight matrix.
#' @param snn.name name of SNN space. If spatial=FALSE and order.cells = NULL, default snn.name will set to 'RNA_snn'. Use SNN space to calculate the cell-cell weight martix.
#' @param order.cells For linear trajetory, input ordered cell names to calculate the cell-cell distance weight matrix. Conflict with sptaial=TRUE and snn.name != NULL.
#' @param cells Cells used for calculate weight matrix. Used with snn graph. In default will use all cells.
#' @param min.cells If a feature can be detect in few than min.cells, will skip to save time. Default is 10.
#' @param features List of features to test. Default is all features with that coverage >= min.cells.
#' @param weight.matrix.name Weight graph name in Seurat object. After this function, the graph can be visited by obj[[weight.matrix.name]]. Default name is "WeightMatrix", if you change the default name, you should specific the new name in \code{\link{{RunBlockCorr}}}.
#' @param prefix Prefix for score and p value names. Default prefix is "moransi". If you change the default name, you should specific the new name in \code{{\link{SetAutoCorrFeatures}}}.
#' @param threads Threads.
#' @param verbose Print log message. Default is TRUE.
#' @param ... parameters pass to \code{\link{GetWeightFromSpatial}}, so it only works if spatial is TRUE.
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
                        min.cells = 10,
                        features = NULL,
                        weight.matrix.name = "WeightMatrix",
                        prefix="moransi",
                        threads=0,
                        verbose = TRUE,
                        ...)
{
  check.par <- 0
  if (!is.null(snn.name)) check.par <- check.par + 1
  if (!is.null(order.cells)) check.par <- check.par + 1
  if (isTRUE(spatial)) check.par <- check.par + 1

  if (check.par > 1) {
    stop("Only support to set one of spatial, snn.name, and order.cells.")
  }
  
  if (check.par == 0) {
    snns <- grep("_snn$", names(object), value=TRUE)
    if (length(snns) == 0) stop("No SNN graph found at object, try to specify snn.name first or run Seurat::FindNeighbors.")
    snn.name <- snns[1]
    if (isTRUE(verbose)) {
      message("Using ", snn.name, " to construct cell weight matrix.")
    }
  }

  if (packageVersion("Seurat") >= numeric_version(as.character(5))) {
    layer <- Layers(object = object, search = layer)
    if (is.null(layer)) {
      abort("No layer found. Please run NormalizeData or RunTFIDF and retry..")
    }
  }
  
  tt <- Sys.time()
  
  assay <- assay %||% DefaultAssay(object = object)
  if (isTRUE(verbose)) {
    message(paste0("Working on assay : ", assay))
  }

  cells <- order.cells %||% cells
  cells <- cells %||% colnames(object)
  cells <- intersect(cells,colnames(object))
  
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
  
  x0 <- GetAssayData1(object, assay = assay, layer = layer)[,cells]
  rs <- rowSums(x0>0)
  features0 <- rownames(x0)[which(rs >= min.cells)]
  features <- intersect(features0, features)
  x0 <- x0[features,]
  x0 <- as(x0, "CsparseMatrix")

  W <- W[cells, cells]
  if (isTRUE(verbose)) {
    message(paste0("Run autocorrelation test for ", length(features), " features."))
  }
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
  if (isTRUE(verbose)) {
    message(paste0("Runtime : ",format(tt)));
  }
  object
}

#' @name SetAutoCorrFeatures
#' @title Set autocorrection by rank and/or score
#' @param object Seurat object
#' @param moransi.min Minimal score for Morans I. In default is 0.
#' @param assay Working assay. If not set, use DefaultAssay(object).
#' @param p.thresh Threshold for p value. Default is 1e-2.
#' @param prefix Prefix name for Moran's index and p value generated by \code{\link{RunAutoCorr}}. Default is "moransi".
#' @param verbose Print log message.
#' @export
SetAutoCorrFeatures <- function(object = NULL,
                                moransi.min = 0,
                                assay = DefaultAssay(object),
                                p.thresh = 0.01,
                                prefix = "moransi",
                                verbose = TRUE)
{
  object0 <- object[[assay]]

  cn <- colnames(object0[[]])
  prefix.p <- paste0(prefix,".pval")
  
  if (prefix %ni% cn | prefix.p %ni% cn) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }

  idx <- which(object0[[prefix]] > moransi.min & object0[[prefix.p]] <= p.thresh)

  if (isTRUE(verbose)) {
    message(paste0(length(idx), " autocorrelated features."))
  }
  object0[["autocorr.variable"]] <- FALSE

  all <- object0[["autocorr.variable"]][[1]]
  all[idx] <- TRUE
  object0[["autocorr.variable"]] <- all
  object[[assay]] <- object0
  
  object
}
#' @name AutoCorrFeatures
#' @title return spatial autocorrlated features
#' @param obejct Seurat object
#' @param assay Working assay.
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
