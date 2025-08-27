#' @export
spatialDistTest <- function(coord = NULL, n = 8) {
  coord <- as.matrix(coord)
  if (nrow(coord) > 10000) {
    coord <- coord[1:10000,]
  }
  m <- .Call("min_Dist_N_Nearest_Neighbors", coord, n)
  m
}

#' @title buildKNN
#' @description Compute KNN with RcppAnnoy.
#' @param k.param K-nearest neighbors, for calculating weight matrix with emb. Default is 20.
#' @importFrom Matrix sparseMatrix
#' @importFrom RANN nn2
#' @export
buildKNN <- function(data, query = data, k.param = 20) {
  nn <- nn2(data, query, k = k.param)
  names(nn) <- c("idx", "dist")
  rownames(nn$idx) <- rownames(query)
  rownames(nn$dist) <- rownames(query)
  nn
}
#' @title buildSNN
#' @description Convert a knn list into a shared nearest neighbors adjacency matrix.
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Default is \code{0}.
#' @export
buildSNN <- function(knn, prune.SNN = 0, diag.value = 1) {
  idx <- knn$idx
  k <- ncol(idx)
  nr <- nrow(idx)
  m <- sparseMatrix(p = c(0,c(1:nr)*k), j = as.vector(t(idx)), x=1)
  m2 <- m %*% t(m)
  m2@x <- m2@x /(k*2-m2@x)
  diag(m2) <- diag.value
  m2@x[m2@x < prune.SNN] <- 0
  m2 <- drop0(m2)
  colnames(m2) <- rownames(knn$idx)
  rownames(m2) <- rownames(knn$idx)
  m2
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
    cells <- rownames(snn)
    snn@x[snn@x <= prune.SNN] <- 0
    snn <- drop0(snn)
    snn <- .Call("norm_by_row", snn)
    snn <- t(snn)
    rownames(snn) <- cells
    colnames(snn) <- cells
    return(snn)
  }

  if (!is.null(pos)) {
    cells <- cells %||% rownames(pos)
    if (is.null(cells)) {
      stop("No rownames for positions.")
    }
    pos <- pos[cells,]
    for (i in 1:ncol(pos)) {
      if (!is.numeric(pos[[i]])) {
        pos[[i]] <- NULL
      }
    }
    pos <- as.matrix(pos)
    knn <- buildKNN(pos, pos, k.param = k.param)
    snn <- buildSNN(knn, prune.SNN = prune.SNN, diag.value = diag.value)
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
  image <- image %||% Images(object)
  if (length(image) == 0) {
    stop("No image found.")
  }

  if (length(image) == 1) {
    emb <- GetTissueCoordinates(object = object, image = image)
    
    if ("cell" %in% colnames(emb)) {
      rownames(emb) <- emb[['cell']]
    }

    W <- GetWeights(pos = emb, diag.value = diag.value, k.param = k.param, prune.SNN=prune.SNN)
    cells <- rownames(emb)
    colnames(W) <- cells
    rownames(W) <- cells
  } else {
    wl <- lapply(image, function(im) {    
      emb <- GetTissueCoordinates(object = object, image = im)
      cells <- emb[['cell']]
      W <- GetWeights(pos = emb, diag.value = diag.value, k.param = k.param, prune.SNN=prune.SNN)
      colnames(W) <- cells
      rownames(W) <- cells
      W
    })
    W <- mergeMatrix(wl)
  }
  return(W)
}
