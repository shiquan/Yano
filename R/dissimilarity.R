#' @importFrom Matrix sparseMatrix
#' @export
RunBlockCorr <- function(object = NULL,
                         bind.name = "gene_name",
                         features = NULL,
                         assay = NULL,
                         bind.assay = NULL,
                         bind.features = NULL,
                         prefix = NULL,
                         #cells = NULL,
                         #order.cells = NULL,
                         feature.types = NULL,                         
                         min.features.per.block = 2,
                         scale.factor = 1e4,
                         weight.matrix.name = "WeightMatrix",
                         #reduction = "pca",
                         mode = 1,
                         #sensitive.mode = FALSE,
                         #spatial = FALSE,
                         #W = NULL,
                         cell.size = NULL,
                         #dims = NULL,
                         #weight.matrix = "RNA_nn",
                         #k.nn = 9,
                         perm=100,
                         threads = 0,
                         block.name = NULL,
                         block.assay = NULL,
                         block.features = NULL
                         )
{
  if (!is.null(block.name)) {
    warning(paste0("'block.name' is deprecated. Set bind.name = '", block.name, "'"))
    bind.name <- block.name
  }

  if (!is.null(block.assay)) {
    warning(paste0("'block.assay' is deprecated. Set bind.assay = '", block.assay, "'"))
    bind.assay <- block.assay
  }

  if (!is.null(block.features)) {
    warning(paste0("'block.features' is deprecated. Set bind.features = '", block.features, "'"))
    bind.features <- block.features
  }

  if (weight.matrix.name %ni% names(object)) {
    stop("No weight matrix found. Perform RunAutoCorr() first.")
  }
  
  tt <- Sys.time()
  cells <- cells %||% colnames(object)
  cells <- intersect(cells, colnames(object))
  
  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))

  prefix <- prefix %||% bind.name
  
  features <- features %||% AutoCorrFeatures(object)
  features <- intersect(features, rownames(object))

  message(paste0("Working on ", length(features), " features."))

  threads <- getCores(threads)

  W <- object[[weight.matrix.name]]
  W <- W[cells, cells]
  tab <- object[[assay]]@meta.features

  if (bind.name %ni% colnames(tab)) {
    stop(paste0("No bind.name found in the feature table of assay ", assay, ". Run LoadEPTanno or LoadVARanno first."))
  }
  
  tab <- tab[tab[[bind.name]] != ".",] # skip unannotated records
  if (!is.null(feature.types) & "type" %in% colnames(tab)) {
    tab <- subset(tab, type %in% feature.types)
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
  
  message(paste0("Processing ", length(blocks), " blocks.."))
  tab <- subset(tab, tab[[bind.name]] %in% blocks)
  
  bind.assay <- bind.assay %||% "tmp.assay"

  x <- GetAssayData(object, assay = assay, slot = "counts")
  
  # cell sizes
  cs <- cell.size %||% colSums(x)  
  #cells <- intersect(cells,names(which(cs > 0)))
  #W <- W[cells, cells]

  #cells <- colnames(W)  
  ncell <- length(cells)

  x <- x[,cells]
  cs <- cs[cells]
  
  ## all.features <- rownames(tab)
  if (bind.assay %ni% names(object)) {
    message("Aggregate counts..")
    x0 <- x[rownames(tab),cells]
    x0 <- as(x0, "TsparseMatrix")
    
    # Aggregate features in the same block
    y <- sparseMatrix(i = match(tab[[bind.name]][x0@i+1], blocks),
                      j = x0@j+1,
                      x = x0@x, dims=c(length(blocks), length(cells)))

    rm(x0)
    rownames(y) <- blocks
    colnames(y) <- cells
  } else {
    message(paste0("Trying to retrieve data from assay ", bind.assay,".."))
    old.assay <- DefaultAssay(object)
    DefaultAssay(object) <- bind.assay
    blocks <- intersect(blocks, rownames(object))
    tab <- subset(tab, tab[[bind.name]] %in% blocks)

    y <- GetAssayData(object, assay = bind.assay, slot = "counts")#[blocks,cells]

    DefaultAssay(object) <- old.assay
    y <- y[,cells]
  }

  features <- intersect(features, rownames(tab))
  tab <- tab[features,]
  bidx <- match(tab[[bind.name]],rownames(y))
  idx <- match(features, rownames(x))
  #cidx <- match(cells, colnames(x))
  
  message(paste0("Test dissimlarity of binding features with ", threads, " threads."))
  gc()
  ta <- .Call("D_test", x, y, W, perm, threads, idx, bidx, cs, scale.factor, mode);
  if (length(ta) == 1) stop(ta[[1]])

  Lx <- ta[[1]]
  Ly <- ta[[2]]
  r <- ta[[3]]
  e <- ta[[4]]
  tval <- ta[[5]]

  names(Lx) <- features
  names(Ly) <- features
  names(r) <- features
  names(e) <- features
  
  pval <- pt(tval, df = perm - 1, lower.tail = FALSE)
  names(pval) <- features
  #padj <- p.adjust(pval, "BH")
  #names(padj) <- features
  tab <- object[[assay]]@meta.features
  tab[[paste0(prefix, ".D")]] <- e[rownames(object)]
  tab[[paste0(prefix, ".r")]] <- r[rownames(object)]
  tab[[paste0(prefix, ".pval")]] <- pval[rownames(object)]
  #tab[[paste0(prefix, ".padj")]] <- padj[rownames(object)]
  object[[assay]]@meta.features <- tab

  rm(ta)
  gc()

  tt <- Sys.time()-tt
  
  message(paste0("Runtime : ",format(tt)));
  object
}
