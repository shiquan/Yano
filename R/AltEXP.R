CheckBindName <- function(object,
                          bind.name,
                          assay = NULL)
{
  old.assay <- DefaultAssay(object)
  assay <- assay %||% old.assay
  DefaultAssay(object) <- assay
  
  meta <- object[[assay]][[]]

  if (bind.name %ni% colnames(meta)) {
    stop("No bind.name found at meta table.")
  }

  sel <- meta[[bind.name]]  
  names(sel) <- rownames(object)

  DefaultAssay(object) <- old.assay
  return(sel)
}

# This function copied from Seurta/R/differential_expression.R
# FindMarkers helper function for cell grouping error checking
ValidateCellGroups <- function(
  object,
  cells.1,
  cells.2,
  min.cells.group = 0
) {
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
}
#
#'@export
PermTest <- function(x, y, cells.1, cells.2, rst, perm = 100, seed = 999, mode = 3, threads = 1, ...)
{
  x <- as(x, "dgCMatrix")
  y <- as(y, "dgCMatrix")
  cells <- colnames(x)
  idx.1 <- match(cells.1, cells)
  idx.2 <- match(cells.2, cells)

  y <- y[rst$bind.feature, ]
  
  df <- .Call("alt_exp", x, y, idx.1, idx.2, mode, perm, threads, seed)

  if (length(df) == 1) {
    stop(df)
  }
  rst$delta <- df[[1]]
  rst$tval <- df[[2]]
  rst$mean <- df[[3]]
  rst$var <- df[[4]]

  rst$pval <- pt(rst$tval, df = perm - 1, lower.tail = FALSE)
  rst$padj <- p.adjust(rst$pval, method = "BH")
  rst
}
DEXSeqTest <- function(x, y, cells.1, cells.2, pesudo.group, rst, mode = 1, threads = 1, ...)
{
  if (!PackageCheck('DEXSeq', error = FALSE)) {
    stop("Install DEXSeq package first. BiocManager::install('DEXSeq')")
  }

  n1 <- length(cells.1)
  n2 <- length(cells.2)

  new.group1 <- rep(1:pesudo.group,(n1/pesudo.group+1), length.out = n1)
  new.group2 <- rep((pesudo.group+1):(pesudo.group*2),(n2/pesudo.group+1), length.out = n2)
  
  new.group <- c(new.group1,new.group2)

  features <- rownames(x)
  x <- x[, c(cells.1, cells.2)]
  x <- as(x,"TsparseMatrix")
  x <- Matrix::sparseMatrix(i=(x@i+1),j=new.group[x@j+1], x=x@x)
  colnames(x) <- paste0("group_",1:(pesudo.group*2))
  rownames(x)<- features

  bind.features <- rownames(y)
  y <- y[, c(cells.1, cells.2)]
  y <- as(y, "TsparseMatrix")
  y <- Matrix::sparseMatrix(i=(y@i+1),j=new.group[y@j+1], x=y@x)
  colnames(y) <- paste0("group_",1:(pesudo.group*2))
  rownames(y) <- bind.features

  rst <- subset(rst, bind.feature %in% bind.features)
  features <- rownames(rst)
  
  meta.tab <- data.frame(row.names=colnames(x),
                         condition=c(rep("cluster1", pesudo.group),
                                     rep("cluster2", pesudo.group)))
  
  if (length(x = features) == 0) {
    warning("No feature pass min.pct threshold; returning NULL.")
    return(NULL)
  }

  x <- as.matrix(x[features,])
  y <- as.matrix(y[rst$bind.feature,])

  stopifnot(identical(dim(x), dim(y)))
  
  if (mode == 2) {
    y <- y - cnt
    y[y<0] <- 0
  }

  if (mode == 3) {
    y <- y + cnt
  }

  meta.tab$condition <- factor(meta.tab$condition)

  feature.names <- paste0("feature_",1:nrow(x))
  names(feature.names) <- features
  rownames(x) <- feature.names
  rownames(y) <- feature.names
  
  dxd <- DEXSeq::DEXSeqDataSet(x,
                               sampleData=meta.tab,
                               groupID = rst$bind.feature,
                               featureID = feature.names,
                               alternativeCountData = y,
                               design= ~sample+exon+condition:exon)

  rm(x)
  rm(y)

  threads <- getCores(threads)
  BPPARAM <- BiocParallel::MulticoreParam(workers = threads)

  results <- DEXSeq::DEXSeq(dxd, BPPARAM=BPPARAM, ...)

  rst$log2fc <- results$log2fold_cluster2_cluster1
  rst$pval <- results$pvalue
  rst$padj <- results$padj

  return(rst)
}
#'
#' @importFrom SeuratObject PackageCheck
#'@export
AlternativeExpressionTest <- function(object,
                                      cells.1 = NULL, cells.2 = NULL,
                                      ident.1 = NULL, ident.2 = NULL,
                                      bind.name = NULL, bind.assay = NULL,
                                      assay = NULL,
                                      layer = "counts",
                                      test.use = c("DEXSeq", "PSI"),
                                      features = NULL,
                                      bind.features = NULL,
                                      min.pct = 0.05,
                                      min.pct.bind.feature = 0.05,
                                      min.cells.group = 3,
                                      pesudo.group = 3,
                                      mode = 1,
                                      threads = 1,
                                      perm = 100,
                                      seed = 999,
                          ...)
{
  if (is.null(bind.name)) {
    stop("No bind.name specified.")
  }
  method <- match.arg(test.use)
  if (method == "PSI") {
    if (is.null(bind.assay)) {
      stop("For PSI test, bind.assay is required.")
    }
    if (mode != 3) {
      message("For PSI test, change mode to 3")
      mode = 3
    }    
  }

  bind.features0 <- CheckBindName(object, bind.name, assay = assay)

  if (!is.null(bind.features)) {
    bind.features <- intersect(bind.features, bind.features0)
  } else {
    bind.features <- bind.features0
  }
  
  if (!is.null(ident.1)) {
    cells.1 <- colnames(object)[which(Idents(object) == ident.1)]
    if (is.null(cells.1)) stop("No ident.1 found")
  }

  if (!is.null(ident.2)) {
    cells.2 <- colnames(object)[which(Idents(object) == ident.2)]
    if (is.null(cells.2)) stop("No ident.2 found")
  }

  if (is.null(cells.1)) stop("No cells.1 specified.")
  
  if (is.null(cells.2)) {
    cells.2 <- setdiff(colnames(object), cells.1)
  }
  
  ValidateCellGroups(object, cells.1, cells.2, min.cells.group)

  features <- features %||% rownames(object)
  
  # filter features
  dat <- GetAssayData1(object, layer=layer)
  x <- dat[features, c(cells.1, cells.2)]
  d1 <- x[, cells.1]
  d2 <- x[, cells.2]
  pct1 <- rowSums(d1>0)/length(cells.1)
  pct2 <- rowSums(d2>0)/length(cells.2)
  pct <- rowSums(x)/length(ncol(x))
  
  rm(d1)
  rm(d2)
  
  rst <- data.frame(feature=features, pct.1 = pct1, pct.2 = pct2)
  rownames(rst) <- features  
  bind.features <- bind.features[features]
  rst$bind.feature <- bind.features
  
  rst <- rst[which(pct >= min.pct & !is.na(rst$bind.feature)),]
  features <- rownames(rst)
  x <- x[features, ]
  
  if (is.null(bind.assay)) {
    x0 <- dat[features,]
    x0 <- as(x0, "TsparseMatrix")

    blocks <- unique(bind.features)
    
    # Aggregate features in the same block
    y <- sparseMatrix(i = match(bind.features[x0@i+1], blocks),
                      j = x0@j+1,
                      x = x0@x, dims=c(length(blocks), ncol(dat)))

    rm(x0)
    colnames(y) <- colnames(dat)
    rownames(y) <- blocks
    y <- y[, c(cells.1, cells.2)]
  } else {
    y <- GetAssayData1(object,assay = bind.assay, layer = layer)
    bind.features0 <- intersect(unique(bind.features), rownames(y))
    b1 <- y[bind.features0, cells.1]
    b2 <- y[bind.features0, cells.2]

    bind.pct1 <- rowSums(b1>0)/length(cells.1)
    bind.pct2 <- rowSums(b2>0)/length(cells.2)
    rm(b1)
    rm(b2)
    
    names(bind.pct1) <- bind.features0
    names(bind.pct2) <- bind.features0
    
    bind.features0 <- bind.features0[which(bind.pct1 >= min.pct.bind.feature & bind.pct2 >= min.pct.bind.feature)]
    bind.features <- bind.features[which(bind.features %in% bind.features0)]

    y <- y[bind.features0, c(cells.1, cells.2)]
    
    features <- names(bind.features)
    rst <- rst[features,]
    rst$bind.feature <- bind.features
    rst$bind.feature.pct.1 <- bind.pct1[bind.features]
    rst$bind.feature.pct.2 <- bind.pct2[bind.features]
  }

  features <- rownames(rst)
  x <- x[features,]
  
  if (method == "DEXSeq") {
    rst <- DEXSeqTest(x, y, cells.1, cells.2, pesudo.group, rst, mode = mode, threads = threads, ...)
  } else if (method == "PSI") {
    rst <- PermTest(x, y, cells.1, cells.2, rst, perm = perm, seed = seed, mode = mode, threads = threads, ...)
  }
  rst
}

#'@importFrom SeuratObject PackageCheck
#' @export
FindAllAEMarkers <- function(object,
                             bind.name = NULL, bind.assay = NULL,
                             layer = "counts",
                             features = NULL,
                             node = NULL,
                             min.pct = 0.1,
                             min.cells.group = 3,
                             min.pct.bind.feature = 0.1,
                             pesudo.group = 3,
                             return.thresh = 1e-2,
                             mode = 1,
                             test.use = c("DEXSeq", "PSI"),
                             threads = 1,
                             perm = 100,
                             seed = 999,
                             ...)

{
  if (is.null(bind.name)) {
    stop("No bind.name specified.")
  }

  ## node paramenter inhert from Seurat v5.0.0, check Seurat::FindAllMarkers and orignal paper for details
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
    if (!PackageCheck('ape', error = FALSE)) {
      stop("Install ape package")
    }
    tree <- Seurat:::Tool(object = object, slot = 'BuildClusterTree')
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- Seurat:::DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
    descendants <- MapVals(
      vec = descendants,
      from = all.children,
      to = tree$tip.label
    )
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(
      node,
      as.numeric(x = setdiff(x = descendants, y = keep.children))
    )
    tree <- ape::drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }

  genes.de <- list()
  messages <- list()
  for (i in 1:length(x = idents.all)) {
    message("Calculating cluster ", idents.all[i])
    genes.de[[i]] <- tryCatch(
      expr = {
        AlternativeExpressionTest(
          object = object,
          ident.1 = if (is.null(x = node)) {
            idents.all[i]
          } else {
            tree
          },
          ident.2 = if (is.null(x = node)) {
            NULL
          } else {
            idents.all[i]
          },
          bind.name = bind.name,
          bind.assay = bind.assay,
          features = features,
          layer = layer,
          min.pct = min.pct,
          min.cells.feature = min.cells.feature,
          min.pct.bind.feature = min.pct.bind.feature,          
          mode = mode,          
          test.use = test.use,
          threads = threads,
          pseudo.group = pseudo.group,
          perm = perm,
          seed = seed,
          ...
        )
      },
      error = function(cond) {
        return(cond$message)
      }
    )
    if (is.character(x = genes.de[[i]])) {
      messages[[i]] <- genes.de[[i]]
      genes.de[[i]] <- NULL
    }
  }
  
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      gde <- subset(x = gde, subset = pval < return.thresh)
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)        
      }
    }
  }
  
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$feature))
  return(gde.all)
}
  
