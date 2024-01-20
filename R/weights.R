## #' Build cell-cell weight matrix on cell lineage.
## #' @param order.cells Predefined cell ranks, used for cell lineage analysis.
## #' @param k.nn k-nearest neighbors.
## #' @param self.weight Diagnoal value in the weight matrix.
## #' @param scale Scale the weight values by row.
## #' @returns A sprase weight matrix.
## #' @importFrom Matrix sparseMatrix
## #' @export
## WeightLineage <- function(order.cells = NULL, k.nn = 9, self.weight = 0, scale = FALSE) {
##   if (is.null(order.cells)) stop("No order.cells specified.")
  
##   n <- length(order.cells)
##   i <- .Call("build_weight_i", n, k.nn)

##   W <- sparseMatrix(i = i + 1, j = rep(1:n, each =k.nn), x=1, dims =c(n,n))

##   W@x[W@x > 1] <- 1
##   diag(W) <- self.weight

##   if (scale) W <- W/rowSums(W)

##   colnames(W) <- order.cells
##   rownames(W) <- order.cells
##   W
## }

#' Calcualte cell-cell weight matrix by one of shared nearest neighbour matrix, spatial locations, cell embedding and linear trajectory.
#' @param snn Shared nearest neighbour matrix, usually can found at object[["RNA_snn"]]
#' @param spatial Tissue coordinates.
#' @param order.cells Predefined cell ranks, used for cell lineage analysis.
#' @param emb Cell dimesional space (PCA/ICA/harmony) 
#' @param k.nn k-nearest neighbors, for emb only.
#' @param diag.value Diagnoal value in the weight matrix.
#' @returns A sparse weight matrix.
#' @importFrom Matrix rowSums
#' @export
GetWeights <- function(snn = NULL,
                       pos = NULL,
                       order.cells = NULL,
                       emb = NULL,
                       k.nn = 20,
                       prune.distance = 20,
                       prune.snn = 1/10,
                       diag.value = 0)
{
  check.par <- 0
  
  if (!is.null(snn)) check.par <- check.par + 1
  if (!is.null(pos)) check.par <- check.par + 1
  if (!is.null(order.cells)) check.par <- check.par + 1
  if (!is.null(emb)) check.par <- check.par + 1
  
  if (check.par != 1) {
    stop("Should only specify one of snn, pos or order.cells.")
  }

  if (!is.null(snn)) {
    diag(snn) <- diag.value
    snn[snn<=prune.snn] <- 0
    snn <- snn/(rowSums(snn)+0.0001)
    return(snn)
  }

  if (!is.null(pos)) {
    pos.dist <- as.matrix(dist(x=pos))
    pos.dist[pos.dist > prune.distance] <- 0
    W <- 1/pos.dist^2
    W[is.na(W)] <- 0
    diag(x = W) <- diag.value
    W <- W/rowSums(W)
    cells <- rownames(pos)

    if (!is.null(cells)) {
      colnames(W) <- cells
      rownames(W) <- cells
    }
    W <- as(W, "CsparseMatrix")
    return(W)
  }

  if (!is.null(emb)) {
    cells <- rownames(emb)
    knn <- buildKNN(emb, k.nn=k.nn)
    snn <- buildSNN(knn)    
    colnames(snn) <- cells
    rownames(snn) <- cells
    W <- GetWeights(snn=snn, diag.value = diag.value, prune.snn = prune.snn)
    return(W)
  }
  
  if (!is.null(order.cells)) {
    pos.dist <- as.matrix(dist(x=c(1:length(order.cells))))
    pos.dist[pos.dist > prune.distance] <- 0
    W <- 1/pos.dist^2
    W[is.na(W)] <- 0
    diag(x = W) <- diag.value
    W <- W/rowSums(W)
    colnames(W) <- order.cells
    rownames(W) <- order.cells
    W <- as(W, "CsparseMatrix")
    return(W)
  }
}

#'
#' @export
GetWeightsFromSNN <- function(object = NULL, name = "RNA_snn", prune.snn = 1/15)
{
  if (name %ni% names(object)) {
    stop(paste0("No ", name, " found at object. Run FindNeighbors on RNA assay first."))
  }

  snn <- object[[name]]
  W <- GetWeights(snn=snn, prune.snn = prune.snn)
  return(W)
}

#' @export
GetWeightsFromSpatial <- function(object = NULL, diag.value = 0, prune.distance = 20, ...) {
  cells <- colnames(object)
  emb <- GetTissueCoordinates(object = object, ...)
  W <- GetWeights(pos = emb, diag.value = diag.value, prune.distance = prune.distance)
  colnames(W) <- cells
  rownames(W) <- cells
  return(W)
}

