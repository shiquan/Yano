#' @title GetWeights
#' @description Calcualte cell-cell weight matrix by one of shared nearest neighbour matrix, spatial locations, cell embedding and linear trajectory.
#' @param snn Shared nearest neighbour graph, usually can found at object[["RNA_snn"]]. This graph can be calculate by Seurat::FindNeighbors().
#' @param pos Tissue coordinates matrix.
#' @param order.cells Predefined cell ranks, used for cell lineage analysis.
#' @param emb Cell dimesional space (PCA/ICA/harmony).
#' @param k.nn K-nearest neighbors, for calculating weight matrix with emb.
#' @param prune.distance Sets the cutoff for cell distance on lineage trajectory (ranked cells) or spatial cooridates (bin distance) when computing the neighborhood overlap for the weight matrix construction. Any edges with values greater than this will be set to 0 and removed from the weight matrix graph. Default is 20, means only calculate weight edges for nearby 20 cells for each cell.
#' @param prune.snn Sets the cutoff for acceptable Jaccard index when computing the neighborhood overlap for the SNN construction. Any edges with values less than or equal to this will be set to 0 and removed from the SNN graph. Essentially sets the stringency of pruning (0 --- no pruning, 1 --- prune everything). Default is 1/15.
#' @param diag.value Diagnoal value in the weight matrix.
#' @param cells Cell list. Default use all cells.
#' @returns A sparse weight matrix.
#' @importFrom Matrix rowSums
#' @export
GetWeights <- function(snn = NULL,
                       pos = NULL,
                       order.cells = NULL,
                       emb = NULL,
                       k.nn = 20,
                       prune.distance = 20,
                       prune.snn = 1/15,
                       diag.value = 0,
                       cells = NULL)
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
    if (!is.null(cells)) {
      snn <- snn[cells, cells]
    }
    snn@x[snn@x <= prune.snn] <- 0
    snn <- snn/(rowSums(snn)+0.0001)
    return(snn)
  }

  if (!is.null(pos)) {
    pos.dist <- as.matrix(dist(x=pos))
    pos.dist[pos.dist > prune.distance] <- 0
    W <- as(pos.dist, "CsparseMatrix")
    W@x <- 1/W@x^2
    diag(x = W) <- diag.value
    W <- W/rowSums(W)
    cells <- rownames(pos)

    if (!is.null(cells)) {
      colnames(W) <- cells
      rownames(W) <- cells
    }
    return(W)
  }

  if (!is.null(emb)) {
    cells <- celles %||% rownames(emb)
    emb <- emb[cells,]
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
    W <- as(pos.dist, "CsparseMatrix")
    W@x <- 1/W@x
    diag(x = W) <- diag.value
    W <- W/rowSums(W)
    colnames(W) <- order.cells
    rownames(W) <- order.cells
    return(W)
  }
}

GetWeightsFromSNN <- function(object = NULL, name = "RNA_snn", prune.snn = 1/15, cells = NULL)
{
  if (name %ni% names(object)) {
    stop(paste0("No ", name, " found at object. Run FindNeighbors on RNA assay first."))
  }

  snn <- object[[name]]
  W <- GetWeights(snn=snn, prune.snn = prune.snn, cells = cells)
  return(W)
}
GetWeightsFromSpatial <- function(object = NULL, diag.value = 0, prune.distance = 20, ...) {
  cells <- colnames(object)
  emb <- GetTissueCoordinates(object = object, ...)
  W <- GetWeights(pos = emb, diag.value = diag.value, prune.distance = prune.distance)
  colnames(W) <- cells
  rownames(W) <- cells
  return(W)
}

