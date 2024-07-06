#'@export
DScore <- function(x = NULL, y = NULL, W = NULL)
{
  d <- .Call("D_score_lite", x, y, W)
  d
}
#' @name RunBlockCorr
#' @title Run spatial dissimilarity test for features and their binding features in parallel.
#' @param object Seurat object
#' @param bind.name Title name for binding features in the meta table. Consider most users start Yano to perform alternative splicing analysis, the default bind.name set to "gene_name".
#' @param features 
#' @param assay
#' @param min.cells
#' @param bind.assay
#' @param bind.features
#' @param min.cells.bind
#' @param prefix
#' @param filter
#' @importFrom Matrix sparseMatrix
#' @export
RunBlockCorr <- function(object = NULL,
                         bind.name = "gene_name",
                         features = NULL,
                         assay = NULL,
                         min.cells = 10,
                         bind.assay = NULL,
                         bind.features = NULL,
                         #cells = NULL,
                         min.cells.bind = 10,
                         prefix = NULL,
                         subset = NULL,
                         min.features.per.block = 1,
                         scale.factor = 1e4,
                         weight.matrix.name = "WeightMatrix",
                         mode = 1,
                         cell.size = NULL,
                         scale = FALSE,
                         perm=100,
                         seed=999,
                         threads = 0,
                         block.name = NULL,
                         block.assay = NULL,
                         block.features = NULL,
                         verbose = TRUE
                         )
{
  if (!is.null(block.name)) {
    if (isTRUE(verbose)) {
      warning(paste0("'block.name' is deprecated. Set bind.name = '", block.name, "'"))
    }
    bind.name <- block.name
  }

  if (!is.null(block.assay)) {
    if (isTRUE(verbose)) {
      warning(paste0("'block.assay' is deprecated. Set bind.assay = '", block.assay, "'"))
    }
    bind.assay <- block.assay
  }

  if (!is.null(block.features)) {
    if (isTRUE(verbose)) {
      warning(paste0("'block.features' is deprecated. Set bind.features = '", block.features, "'"))
    }
    bind.features <- block.features
  }

  if (weight.matrix.name %ni% names(object)) {
    stop("No weight matrix found. Perform RunAutoCorr() first.")
  }
  
  tt <- Sys.time()

  assay <- assay %||% DefaultAssay(object)
  if (isTRUE(verbose)) {
    message(paste0("Working on assay ", assay))
  }

  if (!is.null(bind.assay)) {
    if (bind.assay %ni% names(object)) {
      stop("No bind.assay is found. Make sure you specify the correct assay name.")
    } else {
      if (isTRUE(verbose)) {
        message(paste0("Working on binding assay ", bind.assay))
      }
    }
  }

  prefix <- prefix %||% bind.name
  
  features <- features %||% AutoCorrFeatures(object)
  features <- intersect(features, rownames(object))

  if (isTRUE(verbose)) {
    message(paste0("Working on ", length(features), " features."))
  }
  threads <- getCores(threads)

  W <- object[[weight.matrix.name]]

  #cells <- cells %||% colnames(object)
  cells <- names(which(rowSums(W) > 0))
  ncell <- length(cells)
  
  object0 <- object[[assay]]
  tab <- object0[[]]

  if (bind.name %ni% colnames(tab)) {
    stop(paste0("No bind.name found in the feature table of assay ", assay, ". Run LoadEPTanno or LoadVARanno first."))
  }

  # skip unannotated records
  tab <- tab[tab[[bind.name]] != "." & !is.na(tab[[bind.name]]),] 
  
  ## if (!is.null(feature.types) & "type" %in% colnames(tab)) {
  ##   tab <- subset(tab, type %in% feature.types)
  ## }

  if (!is.null(subset)) {
    tab <- base::subset(tab, subset = subset)
  }
  
  blocks <- names(which(table(tab[[bind.name]]) >= min.features.per.block))

  bind.features <- bind.features %||% blocks
  blocks <- intersect(bind.features, blocks)
  features <- intersect(features, rownames(tab))

  if (length(features) == 0) {
    stop("No features found.")
  }
  
  idx <- match(features, rownames(tab))
  tab0 <- tab[idx,]
  blocks <- intersect(unique(tab0[[bind.name]]), blocks)

  if (isTRUE(verbose)) {
    message(paste0("Processing ", length(blocks), " blocks.."))
  }
  tab <- subset(tab, tab[[bind.name]] %in% blocks)
  features <- intersect(features, rownames(tab))

  bind.assay <- bind.assay %||% "tmp.assay"

  if (bind.assay %ni% names(object)) {

    if (min.features.per.block == 1) {
      if (isTRUE(verbose)) {
        message("No bind.assay specified, update min.features.per.block to 2.")
      }
      min.features.per.block <- 2
      blocks <- names(which(table(tab[[bind.name]]) >= min.features.per.block))
      tab <- subset(tab, tab[[bind.name]] %in% blocks)
      features <- intersect(features, rownames(tab))
    }

    x <- GetAssayData1(object, assay = assay, layer = "counts")
    cs <- cell.size %||% colSums(x)
    
    x <- x[,cells]

    rs <- Matrix::rowSums(x>0)
    idx <- which(rs >= min.cells)
    features1 <- rownames(object)[idx]
    features <- intersect(features, features1)
    
    cs <- cs[cells]
    W <- W[cells,cells]
    if (isTRUE(verbose)) {
      message("Aggregate counts..")
    }
    x0 <- x[rownames(tab),cells]
    x0 <- as(x0, "TsparseMatrix")
    
    # Aggregate features in the same block
    y <- sparseMatrix(i = match(tab[[bind.name]][x0@i+1], blocks),
                      j = x0@j+1,
                      x = x0@x, dims=c(length(blocks), length(cells)))

    rm(x0)
    rownames(y) <- blocks
    colnames(y) <- cells

    norm <- TRUE
  } else {
    if (packageVersion("Seurat") >= numeric_version(as.character(5))) {
      layer <- Layers(object = object0, search = "data")
      if (is.null(layer)) {
        abort("No layer found. Please run NormalizeData or RunTFIDF and retry..")
      }
    }

    x <- GetAssayData1(object, assay = assay, layer = "data")
    x <- x[,cells]
    W <- W[cells,cells]

    rs <- Matrix::rowSums(x>0)
    idx <- which(rs >= min.cells)
    features1 <- rownames(object)[idx]
    features <- intersect(features, features1)
    #tab <- tab[features,]
    blocks <- unique(tab[[bind.name]])
    if (isTRUE(verbose)) {
      message(paste0("Trying to retrieve data from assay ", bind.assay,".."))
    }
    old.assay <- DefaultAssay(object)
    DefaultAssay(object) <- bind.assay
    blocks <- intersect(blocks, rownames(object))
    if (packageVersion("Seurat") >= numeric_version(as.character(5))) {
      layer <- Layers(object = object[[bind.assay]], search = "data")
      if (is.null(layer)) {
        abort(paste0("No layer found. Please run NormalizeData or RunTFIDF for assay ", assay, " and retry.."))
      }
    }

    y <- GetAssayData1(object, assay = bind.assay, slot = "data")
    y <- y[,cells]

    rs <- Matrix::rowSums(y>0)
    idx <- which(rs >= min.cells.bind)
    blocks1 <- rownames(object)[idx]
    blocks <- intersect(blocks, blocks1)

    tab <- subset(tab, tab[[bind.name]] %in% blocks)
    
    cs <- NULL
    norm <- FALSE
    DefaultAssay(object) <- old.assay
  }

  features <- intersect(features, rownames(tab))
  tab <- tab[features,]
  bidx <- match(tab[[bind.name]],rownames(y))
  idx <- match(features, rownames(x))
  if (isTRUE(verbose)) {
    message(paste0("Test dissimlarity of binding features with ", threads, " threads."))
  }
  gc()
  ta <- .Call("D_test", x, y, W, perm, threads, idx, bidx, cs, scale.factor, mode, scale, norm, seed);
  if (length(ta) == 1) stop(ta[[1]])

  Lx <- ta[[1]]
  #Ly <- ta[[2]]
  r <- ta[[3]]
  e <- ta[[4]]
  tval <- ta[[5]]
  mval <- ta[[6]]
  vval <- ta[[7]]
  
  names(Lx) <- features
  #names(Ly) <- features
  names(r) <- features
  names(e) <- features
  names(mval) <- features
  names(vval) <- features
  
  pval <- pt(tval, df = perm - 1, lower.tail = FALSE)
  names(pval) <- features
  tab <- object[[assay]][[]]
  tab[[paste0(prefix, ".D")]] <- e[rownames(object)]
  tab[[paste0(prefix, ".r")]] <- r[rownames(object)]
  tab[[paste0(prefix, ".pval")]] <- pval[rownames(object)]
  tab[[paste0(prefix, ".mean")]] <- mval[rownames(object)]
  tab[[paste0(prefix, ".var")]] <- vval[rownames(object)]

  tab[[paste0(prefix, ".padj")]] <- p.adjust(pval[rownames(object)], method = "BH")
  object0[[colnames(tab)]] <- tab

  object[[assay]] <- object0

  rm(ta)
  gc()

  tt <- Sys.time()-tt
  if (isTRUE(verbose)) {
    message(paste0("Runtime : ",format(tt)));
  }
  object
}

#' @export
cor_dist <- function(x = NULL, y = NULL, W = NULL, perm = 1000, thread = 1)
{
  ta <- .Call("D_distribution_test", x, y, W, perm, thread)
  ta
}
