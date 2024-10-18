#' @title RunAutoCorr
#' @description Calculate spatial autocorrelation (Moran's I) for features in parallel.
#' @param object Seurat object
#' @param assay Working assay
#' @param layer Input data layer, usually be 'data'.
#' @param reduction Cell space used to calculate SNN graph, default is 'pca'.
#' @param dims Dimensions of reduction used to calculate SNN graph.
#' @param prune.SNN This parameter sets the cutoff for the acceptable Jaccard index when computing neighborhood overlap during SNN (Shared Nearest Neighbor) construction. It is passed to Seurat::FindNeighbors. Any edges with Jaccard index values less than or equal to this cutoff will be set to 0 and removed from the SNN graph, effectively controlling the stringency of pruning (with 0 meaning no pruning and 1 meaning everything is pruned). The default value is 1/50, which differs from Seurat's default setting. This is because Seurat’s default is OK for cell clustering but may cause the loss of many sparse features in large cell populations during spatial dissimilarity test. Setting the cutoff to a smaller value can capture more features, but it will also increase computational time.
#' @param k.param Defines k for K-nearest neighbor algorithm
#' @param nn.method nn.method passed to Seurat::FindNeighbors, default is "euclidean".
#' @param n.trees n.trees passed to Seurat::FindNeighbors, default is 50.
#' @param annoy.metric annoy.metric passed to Seurat::FindNeighbors, default is "annoy".
#' @param nn.eps nn.eps passed to Seurat::FindNeighbors, default is 0
#' @param l2.norm L2 normalization. Default is FALSE.
#' @param cells Cells used for calculate weight matrix. Used with snn graph. In default will use all cells.
#' @param min.cells If a feature can be detect in few than min.cells, will skip to save time. Default is 10.
#' @param snn.name name of SNN space. If spatial=FALSE and order.cells = NULL, default snn.name will set to 'RNA_snn'. Use SNN space to calculate the cell-cell weight martix.
#' @param spatial Use spatial coordinate instead of SNN space and linear trajectory to calculate the cell-cell weight matrix.
#' @param order.cells For linear trajetory, input ordered cell names to calculate the cell-cell distance weight matrix. Conflict with sptaial=TRUE and snn.name != NULL.
#' @param weight.method Weight method for distance, default 1/dist^2. Also support average, use mean weight value for nearby cells.
#' @param prune.distance Set the cutoff for neighbors for order cells and spatial coordinates. In default, 50 for order cells, 8 for spatial coordinates.
#' @param features List of features to test. Default is all features with that coverage >= min.cells.
#' @param wm.name Weight matrix/graph name in Seurat object. After this function, the graph can be visited by obj[[wm.name]]. Default name is "RNA_wm", if you change the default name, you should specific the new name in RunBlockCorr.
#' @param prefix Prefix for score and p value names. Default prefix is "moransi". If you change the default name, you should specific the new name in SetAutoCorrFeatures.
#' @param threads Threads.
#' @param verbose Print log message. Default is TRUE.
#' @param ... parameters pass to GetWeightFromSpatial, so it only works if spatial is TRUE.
#' @importFrom methods as
#' @export
RunAutoCorr <- function(object = NULL,
                        assay = NULL,
                        layer = "data",
                        reduction = "pca",
                        dims=1:20,
                        k.param = 20,
                        prune.SNN = 1/50,
                        nn.method = "annoy",
                        n.trees = 50,
                        annoy.metric = "euclidean",
                        nn.eps = 0,
                        l2.norm = FALSE,
                        cells = NULL,
                        min.cells = 10,
                        snn.name = NULL,
                        spatial = FALSE,
                        order.cells = NULL,
                        weight.method = c("dist", "average"),
                        prune.distance = -1,
                        features = NULL,
                        wm.name = NULL,
                        prefix="moransi",
                        threads=0,
                        verbose = TRUE,
                        ...)
{
  check.par <- 0
  if (!is.null(snn.name)) {
    if (snn.name %ni% names(object)) {
      stop(paste0("SNN graph not found, ", snn.name))
    }
    check.par <- check.par + 1
  }
  if (!is.null(order.cells)) check.par <- check.par + 1
  if (isTRUE(spatial)) check.par <- check.par + 1

  if (check.par > 1) {
    stop("Only support to set one of spatial, snn.name, and order.cells.")
  }

  cells <- order.cells %||% cells
  cells <- cells %||% colnames(object)
  cells <- intersect(cells,colnames(object))

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
                        nn.method = nn.method,
                        annoy.metric = annoy.metric,
                        n.trees = n.trees,
                        nn.eps = nn.eps,
                        l2.norm = l2.norm, cache.index = FALSE, verbose = FALSE)
      snn <- ng[['snn']]
      W <- GetWeights(snn = snn, prune.SNN = prune.SNN)
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
  
  if (!is.null(snn.name)) {
    W <- GetWeightsFromSNN(object=object, name = snn.name, prune.SNN = prune.SNN, cells = cells)
  }

  if (!is.null(order.cells)) {
    W <- GetWeights(order.cells = order.cells, prune.distance = prune.distance, weight.method = weight.method)
  }

  if (isTRUE(spatial)) {
    W <- GetWeightsFromSpatial(object = object, prune.distance = prune.distance, weight.method = weight.method, ...)
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

  wm.name <- wm.name %||% paste0(reduction, "_wm")

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
  object0[["autocorr.variable"]] <- all
  object[[assay]] <- object0

  object <- LogSeuratCommand(object)
  
  object
}
#' @title AutoCorrFeatures
#' @description return spatial autocorrlated features
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
