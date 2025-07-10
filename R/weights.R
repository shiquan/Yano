#' @export
spatialDistTest <- function(coord = NULL, n = 8) {
  coord <- as.matrix(coord)
  if (nrow(coord) > 10000) {
    coord <- coord[1:10000,]
  }
  m <- .Call("min_Dist_N_Nearest_Neighbors", coord, n)
  m
}

##
## This function edited from BPCells::knn_annoy(), credit to orignal authors
##
#' @title buildKNN
#' @description Compute KNN with RcppAnnoy.
#' @param k.param K-nearest neighbors, for calculating weight matrix with emb. Default is 20.
#' @param n.trees Number of trees during index build time. More trees gives higher accuracy
#' @export
buildKNN <- function(data, query = data, k.param = 20, metric = c("euclidean", "cosine", "manhattan", "hamming"), n.trees = 50) {
  metric <- match.arg(metric)
  annoy <- switch(metric,
    "euclidean" = new(RcppAnnoy::AnnoyEuclidean, ncol(data)),
    "cosine" = new(RcppAnnoy::AnnoyAngular, ncol(data)),
    "manhattan" = new(RcppAnnoy::AnnoyManhattan, ncol(data)),
    "hamming" = new(RcppAnnoy::AnnoyHamming, ncol(data)),
  )
  for (i in seq_len(nrow(data))) {
    annoy$addItem(i - 1, data[i, ])
  }
  annoy$build(n.trees)

  idx <- matrix(nrow = nrow(query), ncol = k.param)
  dist <- matrix(nrow = nrow(query), ncol = k.param)
  rownames(idx) <- rownames(query)
  rownames(dist) <- rownames(query)
  for (i in seq_len(nrow(query))) {
    res <- annoy$getNNsByVectorList(query[i, ], k.param, -1, include_distances = TRUE)
    idx[i, ] <- res$item + 1
    dist[i, ] <- res$dist
  }
  if (metric == "cosine") dist <- 0.5 * (dist * dist)
  list(idx = idx, dist = dist)
}
#' @title buildKNN
#' @description Convert a knn object into a shared nearest neighbors adjacency matrix.
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 --- prune everything). Default is 1/50.
#' @export
buildSNN <- function(knn, prune.SNN = 1/30) {
  snn <- .Call("knn2snn", knn$idx, prune.SNN)
  snn <- snn + t(snn)
  diag(snn) <- 1
  colnames(snn) <- rownames(knn$idx)
  rownames(snn) <- rownames(knn$idx)
  snn
}
#' @title GetWeights
#' @description Calcualte cell-cell weight matrix by one of shared nearest neighbour matrix, spatial locations, cell embedding and linear trajectory. In default, if no snn/pos/order.cells/emb set, the weight matrix will be generated with PCA.
#' @param snn Shared nearest neighbour graph, usually can found at object[["RNA_snn"]]. This graph can be calculate by Seurat::FindNeighbors().
#' @param pos Manually setup cell coordinates, in a matrix format. This matrix requires at least 2 dimensions.
#' @param order.cells Predefined cell ranks, used for cell lineage analysis.
#' @param emb Cell dimesional space (i.e. PCA/ICA/harmony).
#' @param k.param K-nearest neighbors, for calculating weight matrix with emb/pos/order.cells. Default is 20.
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/30.
## @param nn Nearest neighbors for pos, based on cell coordinates. Default is 8 for spatial coordinate, 20 for lineage trajectory,
#' @param diag.value Diagnoal value in the weight matrix.
#' @param cells Cell list. Default use all cells. If set, pre-defined SNN is not appliable.
#' @returns A column-wise normalised sparse matrix.
#' @importFrom Matrix rowSums sparseMatrix drop0
#' @export
GetWeights <- function(snn = NULL,
                       pos = NULL,
                       order.cells = NULL,
                       emb = NULL,
                       k.param = 20,
                       prune.distance = -1,
                       prune.SNN = 1/30,
                       #nn = -1,
                       diag.value = 0,
                       #weight.method = c("dist", "dist2","average", "label"))
                       cells = NULL)
{
  check.par <- 0
  
  if (!is.null(snn)) check.par <- check.par + 1
  if (!is.null(pos)) check.par <- check.par + 1
  if (!is.null(order.cells)) check.par <- check.par + 1
  if (!is.null(emb)) check.par <- check.par + 1
  
  if (check.par != 1) {
    stop("Should only specify one of snn, pos, emb, or order.cells.")
  }

  if (!is.null(snn)) {
    diag(snn) <- diag.value
    if (!is.null(cells)) {
      snn <- snn[cells, cells]
    }
    snn@x[snn@x <= prune.SNN] <- 0
    snn <- drop0(snn)
    #snn <- snn/(rowSums(snn)+0.0001)
    snn <- snn/rowSums(snn)
    snn <- t(snn)
    return(snn)
  }

  if (!is.null(pos)) {
    cells <- cells %||% rownames(pos)
    if (is.null(cells)) {
      stop("No rownames for positions.")
    }
    pos <- pos[cells,]
    knn <- buildKNN(pos, pos, k.param = nn)
    snn <- buildSNN(knn, prune.SNN = prune.SNN)
    colnames(snn) <- cells
    rownames(snn) <- cells
    W <- GetWeights(snn=snn, diag.value = diag.value, prune.SNN = prune.SNN)
    return(W)
  }

  if (!is.null(emb)) {
    cells <- cells %||% rownames(emb)
    emb <- emb[cells,]
    knn <- buildKNN(emb, k.param=k.param)
    snn <- buildSNN(knn, prune.SNN = prune.SNN)    
    colnames(snn) <- cells
    rownames(snn) <- cells
    W <- GetWeights(snn=snn, diag.value = diag.value, prune.SNN = prune.SNN)
    return(W)
  }
  
  if (!is.null(order.cells)) {
    pos.dist <- as.matrix(dist(x=c(1:length(order.cells))))
    pos.dist[pos.dist > k.param] <- 0
    W <- as(pos.dist, "CsparseMatrix")
    diag(x = W) <- diag.value
    W <- W/rowSums(W)
    W[is.na(W)] <- 0
    W <- as(W, "CsparseMatrix")
    W <- drop0(W)
    W <- t(W)
    colnames(W) <- order.cells
    rownames(W) <- order.cells
    return(W)
  }
}

GetWeightsFromSNN <- function(object = NULL, snn = "RNA_snn", prune.SNN = 1/50, cells = NULL)
{
  if (snn %ni% names(object)) {
    stop(paste0("No ", snn, " found at object. Run FindNeighbors on RNA assay first."))
  }
  snn.graph <- object[[snn]]
  W <- GetWeights(snn=snn.graph, prune.SNN = prune.SNN, cells = cells)
  return(W)
}
#'@export
GetWeightsFromSpatial <- function(object = NULL, diag.value = 0, k.param = 20, image = NULL, prune.SNN = 1/30) {
  wl <- lapply(image, function(im) {    
    emb <- GetTissueCoordinates(object = object, image = im)
    W <- GetWeights(pos = emb, diag.value = diag.value, k.param = k.param, prune.SNN=prune.SNN)
    cells <- rownames(emb)
    colnames(W) <- cells
    rownames(W) <- cells
    W
  })
  W <- mergeMatrix(wl)
  return(W)
}
