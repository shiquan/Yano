#' @export
makeMetaCells <- function(object, features = NULL, ncell = 5000L, assay = NULL, leverage.name = "leverage.score", seed = 123)
{
  ncols <- ncol(object)
  
  if (ncell < 1) {
    ncell <- ncell * ncols
  } else if (ncell > ncols) {
    ncell <- ncols
  }
  
  assay <- assay %||% DefaultAssay(object)
  ## for v5
  assay <- match.arg(arg = assay, choices = Assays(object = object))
  features <- features %||% VariableFeatures(object)
  object <- LeverageScore(object = object, assay = assay, features = features, var.name = leverage.name, nsketch = ncell)
  leverage.score <- object[[leverage.name]]

  ## the following coding edited from Seurat/R/sketch.R
  layers.data <- Layers(object = object[[assay]], search = 'data')

  cells <- lapply(
    X = seq_along(along.with = layers.data),
    FUN = function(i, seed) {
      set.seed(seed = seed)
      lcells <- Cells(x = object[[assay]], layer = layers.data[i])
      if (length(x = lcells) < ncell) {
        return(lcells)
      }
      return(sample(
        x = lcells,
        size = ncell,
        prob = leverage.score[lcells,]
      ))
    },
    seed = seed
  )

  cells <- unlist(cells)

  cells
}
#' @export
addMetaCells <- function(object, assays = NULL, meta.cells = NULL, reduction = "pca", dims = 1:30, meta.name = "meta_cells", k.param = 1, ncell = 5000L, features = NULL, seed = 123)
{
  if (is.null(meta.cells)) {
    #stop("No meta.cells, use makeMetaCells select meta cells.")
    meta.cells <- makeMetaCells(object, ncell = ncell, features = features, seed = seed)
  }

  data <- Reductions(object, reduction)[[]]

  gp <- FindNeighbors(data[meta.cells,dims], query = data[,dims], return.neighbor = TRUE, k.param=k.param, compute.SNN = FALSE)

  object[[meta.name]] <- paste0("metacell-", gp@nn.idx[,1])
  #orig.name <- meta.cells
  #names(orig.name) <- paste0("metacell-", 1:length(orig.name))
  # object[['meta_nn']] <- gp
  
  #object0 <- AggregateExpression(object, groups.by = meta.name, return.seurat = TRUE)
  #object0$orig.name <- orig.name[colnames(object0)]
  object
}
RefinedDist <- function(dist, coord, dd)
{
  .Call("refinedDist", dist, coord, dd)
}

minIdx <- function(dist, invert = FALSE) {
  .Call("minIdx", dist, invert)
}
#' @export
addSpatialMetaCells <- function(object, assays = NULL, meta.cells = NULL, reduction = "pca", dims =1:50, k.param = 30, meta.name = "meta_cells", slice = NULL, use.reduction = NULL, features = NULL, threads = 0, verbose = TRUE, ncell = 5000L, seed = 123)
{
  threads <- getCores(threads)

  if (is.null(meta.cells)) {
    if (verbose) {
      message("Make meta cells.")
    }
    meta.cells <- makeMetaCells(object, features = features, ncell = ncell, seed = seed)
    #stop("No meta.cells, use makeMetaCells select meta cells.")
  }

  object[['metacell_anchor']] <- FALSE
  object[['metacell_anchor']][meta.cells,] <- TRUE

  coord <- NULL
  if (!is.null(slice) & !is.null(use.reduction)) {
    stop("slice and use.reduction cannot specified at the same.")
  }
  
  if (is.null(slice) & is.null(use.reduction)) {
      slice <- slice %||% Images(object)
      if (is.null(slice)) {
        stop("No slice name found or use.reduction specified .")
      }
  }
  if (!is.null(slice)) {
    coord <- GetTissueCoordinates(object, slice)
    coord <- as.matrix(coord[,c(1,2)])
  } else {
    coord <- Reductions(object, use.reduction)[[]]
    coord <- as.matrix(coord[,c(1,2)])
  }
  
  if (reduction %in% Reductions(object)) {
    data <- Loadings(object[[reduction]])
    data <- data[,dims]
  } else {
    features <- features %||% VariableFeatures(object)
    object0 <- object[,meta.cells]
    object0 <- ScaleData(object0, verbose = verbose)
    object0 <- RunPCA(object0, verbose = verbose)
    data <- ProjectCellEmbeddings(object, object0, dims=dims, verbose = verbose)
    rm(object0)
  }
  gp <- FindNeighbors(data[meta.cells, ], query = data, return.neighbor = TRUE, k.param=k.param, compute.SNN = FALSE)
  # object[["gp"]] <- gp
  ## cells <- rownames(data)
  ## m <- .Call("query_nn", data, match(meta.cells, cells), k.param, threads, dims[2]);
  ## rownames(m) <- cells
  ## colnames(m) <- cells

  dd <- spatialDistTest(coord)

  m <- sparseMatrix(i = rep(1:nrow(gp@nn.idx), each = ncol(gp@nn.idx)),
                    j = as.vector(t(gp@nn.idx)),
                    x = as.vector(t(gp@nn.dist)),
                    dims = c(nrow(gp@nn.idx), nrow(gp@nn.idx)))
  cells <- gp@cell.names
  rownames(m) <- cells
  colnames(m) <- c(meta.cells, setdiff(cells, meta.cells))
  m <- m[cells, cells]
  
  m <- t(m)
  object[["nn"]] <- as(m, "Graph")
  
  coord <- coord[cells,]
  RDD <- RefinedDist(m, coord, dd)
  colnames(RDD) <- cells
  rownames(RDD) <- cells
  object[["spnn"]] <- as(RDD, "Graph")
  rm(m)
  
  idx <- minIdx(RDD, FALSE)
  object[[meta.name]] <- as.vector(paste0("metacell-",idx))
  object
}
#' @export
transferLabels <- function(object, labels, cluster.name = "meta.clusters", meta.name = "meta_cells")
{
  if (is.null(labels)) {
    stop("No labels.")
  }

  object[[cluster.name]] <- as.vector(labels[object[[meta.name]][[1]]])
  object
}

#' @export
transferLabelsFromMeta <- function(object, meta, meta.cluster.name = NULL, cluster.name = "meta.clusters", meta.name = "meta_cells", orig.name = NULL)
{
  if (is.null(meta.cluster.name)) {
    stop("No meta.cluster.name for meta object specified.")
  }

  if (meta.cluster.name %ni% colnames(meta[[]])) {
    stop("No found meta.cluster.name in the meta table.")
  }

  clusters <- meta[[meta.cluster.name]][[1]]
  
  if (is.null(orig.name)) {
    names(clusters) <- colnames(meta)
  } else {
    names(clusters) <- meta[[orig.name]][[1]]
  }

  object <- transferLabels(object, clusters, cluster.name, meta.name)

  object
}
