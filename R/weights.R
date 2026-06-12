#' @title spatialDistTest
#' @description Compute minimum distances to N nearest neighbors from spatial coordinates.
#' @param coord A matrix of spatial coordinates.
#' @param n Number of nearest neighbors. Default is 8.
#' @return A matrix of distances.
#' @export
spatialDistTest <- function(coord = NULL, n = 8) {
  coord <- as.matrix(coord)
  if (n < 1) stop("n must be at least 1.")
  if (n > nrow(coord)) stop("n cannot exceed the number of coordinates.")
  if (nrow(coord) > 10000) {
    warning("More than 10000 coordinates, truncating to first 10000 rows.")
    coord <- coord[1:10000,]
  }
  m <- .Call("min_Dist_N_Nearest_Neighbors", coord, n)
  m
}

##
## This function edited from BPCells::knn_annoy(), credit to original authors
##
#' @title buildKNN
#' @description Compute KNN with RcppAnnoy.
#' @param data Reference data matrix (rows = items, columns = features).
#' @param query Query data matrix. Defaults to data.
#' @param k.param Number of nearest neighbors. Default is 20.
#' @param metric Distance metric. One of "euclidean", "cosine", "manhattan", or "hamming".
#' @param n.trees Number of trees during index build time. More trees gives higher accuracy. Default is 50.
#' @return A list with \code{idx} (neighbor indices) and \code{dist} (neighbor distances).
#' @export
buildKNN <- function(data, query = data, k.param = 20, metric = c("euclidean", "cosine", "manhattan", "hamming"), n.trees = 50) {
  data <- as.matrix(data)
  query <- as.matrix(query)
  if (nrow(data) == 0 || ncol(data) == 0) stop("data must have at least 1 row and 1 column.")
  if (nrow(query) == 0) stop("query must have at least 1 row.")
  if (ncol(data) != ncol(query)) {
    stop("data and query must have the same number of columns.")
  }
  if (nrow(data) <= k.param) {
    k.param <- max(1, nrow(data) - 1)
  }
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

  idx <- matrix(0L, nrow = nrow(query), ncol = k.param)
  dist <- matrix(0, nrow = nrow(query), ncol = k.param)
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
#' @title buildSNN
#' @description Convert a KNN object into a shared nearest neighbors adjacency matrix.
#' @param knn A list with \code{idx} matrix, as returned by \code{buildKNN}.
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 --- prune everything). Default is 1/30.
#' @return A sparse symmetric SNN matrix (\code{dgCMatrix}).
#' @export
buildSNN <- function(knn, prune.SNN = 1/30) {
  if (is.null(knn$idx) || !is.matrix(knn$idx)) {
    stop("knn must be a list with an 'idx' matrix, as returned by buildKNN.")
  }
  idx <- knn$idx
  storage.mode(idx) <- "integer"
  snn <- .Call("knn2snn", idx, prune.SNN)
  snn <- snn + t(snn)
  diag(snn) <- 1
  colnames(snn) <- rownames(knn$idx)
  rownames(snn) <- rownames(knn$idx)
  snn
}
#' @title GetWeights
#' @description Calculate cell-cell weight matrix from one of: shared nearest neighbour graph, spatial positions, cell embeddings, or linear order. Exactly one source must be specified.
#' @param snn Shared nearest neighbour graph, usually found at object[["RNA_snn"]]. This graph can be calculated by Seurat::FindNeighbors().
#' @param pos Manually set cell positions, as a data.frame or matrix with at least 2 numeric columns and rownames as cell names.
#' @param order.cells Predefined cell ranks, used for cell lineage analysis.
#' @param emb Cell embedding space (i.e. PCA/ICA/harmony).
#' @param k.param K-nearest neighbors, for calculating weight matrix with emb/pos/order.cells. Default is 20.
#' @param prune.SNN Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Default is 1/30.
#' @param diag.value Diagonal value in the weight matrix. Default is 0.
#' @param cells Cell list to subset. Default uses all cells.
#' @returns A column-wise normalised sparse matrix.
#' @importFrom Matrix rowSums sparseMatrix drop0
#' @export
GetWeights <- function(snn = NULL,
                       pos = NULL,
                       order.cells = NULL,
                       emb = NULL,
                       k.param = 20,
                       prune.SNN = 1/30,
                       diag.value = 0,
                       cells = NULL)
{
  if (sum(!is.null(snn), !is.null(pos), !is.null(order.cells), !is.null(emb)) != 1) {
    stop("Should only specify one of snn, pos, emb, or order.cells.")
  }

  if (!is.null(snn)) {
    if (!inherits(snn, "sparseMatrix")) {
      stop("snn must be a sparse matrix (dgCMatrix).")
    }
    # make a copy to avoid mutating the caller's matrix
    snn <- as(snn, "CsparseMatrix")
    diag(snn) <- diag.value
    if (!is.null(cells)) {
      cells <- intersect(cells, rownames(snn))
      if (length(cells) == 0) return(NULL)
      snn <- snn[cells, cells]
    }
    snn@x[snn@x <= prune.SNN] <- 0
    snn <- drop0(snn)
    snn <- snn/rowSums(snn)
    snn <- t(snn)
    return(snn)
  }

  if (!is.null(pos)) {
    cells <- cells %||% rownames(pos)
    if (is.null(cells)) {
      stop("No rownames for positions.")
    }
    cells <- intersect(cells, rownames(pos))
    if (length(cells) == 0) {
      return(NULL)
    }
    pos <- pos[cells, , drop = FALSE]
    pos <- pos[, apply(pos, 2, is.numeric), drop = FALSE]
    if (ncol(pos) < 2) stop("pos must have at least 2 numeric columns.")
    pos <- as.matrix(pos)
    knn <- buildKNN(pos, pos, k.param = k.param)
    snn <- buildSNN(knn, prune.SNN = prune.SNN)
    colnames(snn) <- cells
    rownames(snn) <- cells
    W <- GetWeights(snn=snn, diag.value = diag.value, prune.SNN = prune.SNN)
    return(W)
  }

  if (!is.null(emb)) {
    cells <- cells %||% rownames(emb)
    cells <- intersect(cells, rownames(emb))
    if (length(cells) == 0) return(NULL)
    emb <- emb[cells, , drop = FALSE]
    knn <- buildKNN(emb, k.param=k.param)
    snn <- buildSNN(knn, prune.SNN = prune.SNN)    
    colnames(snn) <- cells
    rownames(snn) <- cells
    W <- GetWeights(snn=snn, diag.value = diag.value, prune.SNN = prune.SNN)
    return(W)
  }
  
  if (!is.null(order.cells)) {
    n <- length(order.cells)
    if (n <= 1) stop("order.cells must have at least 2 elements.")
    k.param <- min(k.param, n - 1)
    nz <- n * (2 * k.param + 1) - k.param * (k.param + 1)
    I <- integer(nz); J <- integer(nz); X <- numeric(nz)
    pos <- 1L
    for (d in 0:k.param) {
      if (d > 0) {
        len <- n - d
        idx <- pos:(pos + 2 * len - 1L)
        I[idx] <- c(1:len, (1 + d):n)
        J[idx] <- c((1 + d):n, 1:len)
        X[idx] <- rep(d, 2 * len)
        pos <- pos + 2L * len
      } else {
        idx <- pos:(pos + n - 1L)
        I[idx] <- 1:n
        J[idx] <- 1:n
        X[idx] <- 0
        pos <- pos + n
      }
    }
    stopifnot(pos - 1L == nz)
    W <- sparseMatrix(i = I, j = J, x = X, dims = c(n, n))
    diag(W) <- diag.value
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

#' @title GetWeightsFromSNN
#' @description Calculate cell-cell weight matrix from a precomputed SNN graph stored in the Seurat object.
#' @param object Seurat object.
#' @param snn Name of the SNN graph. Default is "RNA_snn".
#' @param prune.SNN Jaccard index cutoff for SNN pruning. Default is 1/50.
#' @param cells Cell list to subset. Default uses all cells.
#' @return A column-wise normalized sparse weight matrix.
#' @export
GetWeightsFromSNN <- function(object = NULL, snn = "RNA_snn", prune.SNN = 1/50, cells = NULL)
{
  if (snn %ni% names(object)) {
    stop(paste0("No ", snn, " found at object. Run FindNeighbors on RNA assay first."))
  }
  snn.graph <- object[[snn]]
  if (is.null(snn.graph)) stop(paste0("Graph ", snn, " is NULL."))
  W <- GetWeights(snn=snn.graph, prune.SNN = prune.SNN, cells = cells)
  return(W)
}
# Parse coordinate data.frames to a standard format
# Extracts x/y as numeric, sets rownames from 'cell' column or existing rownames
.parseCoords <- function(emb) {
  if (!all(c('x', 'y') %in% colnames(emb))) {
    stop("coords must have 'x' and 'y' columns.")
  }
  if ('cell' %in% colnames(emb)) {
    rn <- as.character(emb[['cell']])
    emb <- data.frame(x = as.numeric(emb[["x"]]), y = as.numeric(emb[["y"]]))
    rownames(emb) <- rn
  } else {
    rn <- rownames(emb)
    if (is.null(rn)) stop("coords must have rownames or a 'cell' column.")
    emb <- data.frame(x = as.numeric(emb[["x"]]), y = as.numeric(emb[["y"]]))
    rownames(emb) <- rn
  }
  emb
}

#' @title GetWeightsFromSpatial
#' @description Calculate cell-cell weight matrix from spatial coordinates.
#' @param object Seurat object.
#' @param diag.value Diagonal value in the weight matrix. Default is 0.
#' @param k.param K-nearest neighbors. Default is 20.
#' @param image Name of the spatial image. If NULL, uses all images in the object.
#' If multiple images are available, each generates a separate weight matrix which
#' are then merged.
#' @param prune.SNN Jaccard index cutoff for SNN pruning. Default is 1/30.
#' @param cells Cell list to subset. Default uses all cells.
#' @return A column-wise normalized sparse weight matrix.
#' @export
GetWeightsFromSpatial <- function(object = NULL, diag.value = 0, k.param = 20, image = NULL, prune.SNN = 1/30, cells = NULL) {
  image <- image %||% Images(object)
  if (length(image) == 0) {
    stop("No image found.")
  }

  if (length(image) == 1) {
    emb <- GetTissueCoordinates(object = object, image = image)
    emb <- .parseCoords(emb)
    W <- GetWeights(pos = emb, diag.value = diag.value, k.param = k.param, prune.SNN = prune.SNN, cells = cells)
  } else {
    wl <- lapply(image, function(im) {
      emb <- GetTissueCoordinates(object = object, image = im)
      emb <- .parseCoords(emb)
      GetWeights(pos = emb, diag.value = diag.value, k.param = k.param, prune.SNN = prune.SNN, cells = cells)
    })
    wl <- Filter(Negate(is.null), wl)
    if (length(wl) == 0) return(NULL)
    W <- mergeMatrix(wl)
  }
  return(W)
}
#' @title GetWeightsFromCoords
#' @description Calculate cell-cell weight matrix from custom spatial coordinates.
#' @param coords A data.frame with columns 'x' and 'y', or a list of such data.frames for multiple slices.
#' The data.frame should have rownames as cell names, or a 'cell' column specifying cell identities.
#' @param diag.value Diagonal value in the weight matrix. Default is 0.
#' @param k.param K-nearest neighbors for building the SNN graph. Default is 30.
#' @param prune.SNN Jaccard index cutoff for SNN pruning. Edges with values below this are removed. Default is 1/30.
#' @param cells Cell list to subset. Default uses all cells in coords.
#' @return A column-wise normalized sparse weight matrix.
#' @export
GetWeightsFromCoords <- function(coords = NULL, diag.value = 0, k.param = 30, prune.SNN = 1/30, cells = NULL) {
  if (is.data.frame(coords)) {
    emb <- .parseCoords(coords)
    W <- GetWeights(pos = emb, diag.value = diag.value, k.param = k.param, prune.SNN = prune.SNN, cells = cells)
  } else if (is.list(coords)) {
    wl <- lapply(coords, function(x) {
      emb <- .parseCoords(x)
      GetWeights(pos = emb, diag.value = diag.value, k.param = k.param, prune.SNN = prune.SNN, cells = cells)
    })
    wl <- Filter(Negate(is.null), wl)
    if (length(wl) == 0) return(NULL)
    W <- mergeMatrix(wl)
  } else {
    stop("coords must be a data.frame or a list of data.frames.")
  }
  return(W)
}
