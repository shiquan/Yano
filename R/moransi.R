#' @title RunAutoCorr
#' @description Calculate spatial autocorrelation (Moran's I) for features in parallel. Autocorrelated features are labeled with SetAutoCorrFeatures() automatically. 
#' @param object Seurat object
#' @param assay Working assay
#' @param layer Input data layer, usually be 'data'.
#' @param snn Name of SNN graph. Use the SNN graph to calculate the cell-cell weight martix. The SNN graph usually be generated by Seurat::FindNeighbors().
#' @param reduction If snn is not set, the cell space will be used to calculate SNN graph. Default is 'pca'.
#' @param dims Dimensions of reduction used to calculate SNN graph.
#' @param prune.SNN This parameter sets the cutoff for the acceptable Jaccard index. It is passed to Seurat::FindNeighbors. Default is 1/30.
#' @param k.param Defines k for K-nearest neighbor algorithm
#' @param cells Cells used for calculate weight matrix. Used with snn graph. In default will use all cells.
#' @param min.cells If a feature can be detect in few than min.cells, will skip to save time. Default is 10.
#' @param image Name of 'SpatialImage' object to get coordinates for. If set, use spatial coordinate to calculate the cell-cell weight matrix.
#' @param order.cells For linear trajetory, input ordered cell names to calculate the cell-cell distance weight matrix. Conflict with sptaial=TRUE and snn.name != NULL.
#' @param weight.method Weight method for distance, default 1/dist^2. Also support average, use mean weight value for nearby cells.
#' @param nn Set the cutoff for nearest neighbors for order cells and spatial coordinates. In default, 50 for order cells, 8 for spatial coordinates.
#' @param features List of features to test. Default is all features with that coverage >= min.cells.
#' @param wm.name Weight matrix/graph name in Seurat object. After this function, the graph can be visited by obj[[wm.name]]. Default name is "$reduction$_wm" for reduction, "trajectory_wm" for order.cells, and "spatial_wm" for spatial coordinate.
#' @param prefix Prefix for score and p value names. Default prefix is "moransi". If you change the default name, you should specific the new name in SetAutoCorrFeatures.
#' @param threads Threads.
#' @param verbose Print log message. Default is TRUE.
#' @param name Title name to label autocorrelation features. Default is 'autocorr.variable'.
#' @param ... parameters pass to FindNeighbors.
#' @importFrom methods as
#' @export
RunAutoCorr <- function(object = NULL,
                        assay = NULL,
                        layer = "data",
                        reduction = "pca",
                        dims=1:20,
                        k.param = 20,
                        prune.SNN = 1/30,
                        ident = NULL,
                        cells = NULL,
                        min.cells = 10,
                        snn = NULL,
                        order.cells = NULL,
                        image = NULL,
                        weight.method = c("dist", "average"),
                        nn = -1,
                        features = NULL,
                        wm.name = NULL,
                        prefix="moransi",
                        threads=0,
                        verbose = TRUE,
                        name = "autocorr.variable",
                        ...)
{
  check.par <- 0
  if (!is.null(snn)) {
    if (snn %ni% names(object)) {
      stop(paste0("SNN graph not found, ", snn))
    }
    check.par <- check.par + 1
  }
  if (!is.null(order.cells)) check.par <- check.par + 1
  if (!is.null(image)) check.par <- check.par + 1

  if (check.par > 1) {
    stop("Only support to set one of image, snn, and order.cells.")
  }

  if (!is.null(ident) && !is.null(cells)) {
    stop("ident is conflict with cells")
  }

  if (!is.null(order.cells) && !is.null(cells)) {
    stop("order.cells is conflict with cells")
  }

  if (!is.null(order.cells) && !is.null(ident)) {
    stop("order.cells is conflict with cells")
  }

  if (!is.null(ident)) {
    cells <- colnames(object)[which(Idents(object) == ident)]
    if (is.null(cells)) stop("No ident found")
  } else {
    cells <- order.cells %||% cells
    cells <- cells %||% colnames(object)
    cells <- intersect(cells,colnames(object))
  }
  
  if (check.par == 0) {
    if (reduction %ni% Reductions(object)) {
      stop(paste0("No reduction ", reduction, " found. You should perform RunPCA first."))
    }
    data.use <- Embeddings(object[[reduction]])
    data.use <- data.use[cells, dims]
    ng <- FindNeighbors(object = data.use,
                        k.param = k.param,
                        compute.SNN = TRUE,
                        prune.SNN = prune.SNN,
                        cache.index = FALSE, verbose = FALSE, ...)
    snn.graph <- ng[['snn']]
    W <- GetWeights(snn = snn.graph, prune.SNN = prune.SNN)
    wm.name <- wm.name %||% paste0(reduction, "_wm")
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
  
  features <- features %||% rownames(object)
  features <- intersect(rownames(object),features)
  
  threads <- getCores(threads)
  
  if (!is.null(snn)) {
    W <- GetWeightsFromSNN(object=object, snn = snn, prune.SNN = prune.SNN, cells = cells)

    wm.name <- wm.name %||% gsub("(.*)_snn", "\\1_wm", snn)
  }

  if (!is.null(order.cells)) {
    W <- GetWeights(order.cells = order.cells, nn = nn, weight.method = weight.method)
    wm.name <- wm.name %||% "trajectory_wm"
  }

  if (!is.null(image)) {
    W <- GetWeightsFromSpatial(object = object, nn = nn, weight.method = weight.method, image = image)
    wm.name <- wm.name %||% "spatial_wm"
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
  
  object[[wm.name]] <- W
  
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

  ## trim logs
  cells <- head(cells)
  features <- head(features)
  
  object <- LogSeuratCommand(object)
  
  tt <- Sys.time()-tt
  if (isTRUE(verbose)) {
    message(paste0("Runtime : ",format(tt)))
  }
  
  object <- SetAutoCorrFeatures(object, name = name, prefix = prefix)
  
  object
}

#' @title SetAutoCorrFeatures
#' @description Set autocorrection by Moran's I and/or p value.
#' @param object Seurat object
#' @param moransi.min Minimal score for Morans I. In default is 0.
#' @param assay Working assay. If not set, use DefaultAssay(object).
#' @param p.thresh Threshold for p value. Default is 1e-2.
#' @param prefix Prefix name for Moran's index and p value generated by RunAutoCorr. Default is "moransi".
#' @param verbose Print log message.
#' @export
SetAutoCorrFeatures <- function(object = NULL,
                                moransi.min = 0,
                                assay = DefaultAssay(object),
                                p.thresh = 0.01,
                                prefix = "moransi",
                                name = "autocorr.variable",
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
  names(all) <- rownames(object0)
  object0[[name]] <- all
  object[[assay]] <- object0

  object <- LogSeuratCommand(object)
  
  object
}
#' @title AutoCorrFeatures
#' @description return spatial autocorrlated features
#' @param obejct Seurat object
#' @param assay Working assay.
#' @export
AutoCorrFeatures <- function(object = NULL, assay = NULL, name = "autocorr.variable")
{
  assay <- assay %||% DefaultAssay(object)

  object0 <- object[[assay]]

  if (name %ni% colnames(object0[[]])) {
    stop("No autocorrelation flag found, run SetAutoCorrFeatures() first.")
  }
  rownames(object)[which(object0[[name]][[1]])]
}
