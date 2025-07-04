DScore <- function(x = NULL, y = NULL, W = NULL)
{
  d <- .Call("D_score_lite", x, y, W)
  d
}

#
#' @importFrom methods is
#
IdentToCells <- function(
  object,
  ident,
  cellnames.use
) {
  #
  if (is.null(x = ident)) {
    stop("Please provide ident")
  } else if ((length(x = ident) == 1 && ident[1] == 'clustertree') || is(object = ident, class2 = 'phylo')) {
    tree <- if (is(object = ident, class2 = 'phylo')) {
      ident
    } else {
      Tool(object = object, slot = 'BuildClusterTree')
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident'")
    }
    ident <- tree$tip.label[GetLeftDescendants(tree = tree, node = ident)]
  }
  if (length(x = as.vector(x = ident)) > 1 &&
        any(as.character(x = ident) %in% cellnames.use)) {
    bad.cells <- cellnames.use[which(x = !as.character(x = ident) %in% cellnames.use)]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident <- WhichCells(object = object, idents = ident)
  }
  return(ident)
}

SDT <- function(x, y, idx, bidx, W, cs, threads, perm, scale.factor, mode, scale, norm, seed, debug)
{
  cells <- colnames(W)
  ta <- .Call("D_test", x[,cells], y[,cells], W, 1, perm, threads, idx, bidx, cs[cells], scale.factor, mode, scale, norm, seed, debug)
  ta
}

#' @title RunBlockCorr
#' @description Run spatial dissimilarity test for features and their binding features in parallel.
#' @param object Seurat object
#' @param bind.name The title of binding features in the meta table is important, as most users begin Yano to perform alternative splicing analysis. By default, the `bind.name` is set to "gene_name".
#' @param features Vector of features to calculate. Default is AutoCorrFeatures(object).
#' @param assay Work assay.
#' @param bind.assay Name of binding assay. If the binding assay is not set, raw counts of features from the same block will be aggregated first, followed by normalization.
#' @param bind.features List of bind features. If set, bind features will be subset with this list.
#' @param min.cells.bind Binding features detected in few than minimal number of cells will be skipped. Default is 10.
#' @param prefix Prefix name for output scores and values. Default is same with bind.name.
## @param subset Rules for subsetting meta table before selecting features to perform test.
#' @param scale.factor The scale factor used to normalize counts is set to a default of 1e4. The raw counts from "counts" layer are further normalized by sample size and scale factor for the spatial dissimilarity test.
#' @param mode Test mode. Mode 1: Compares X (test feature) vs Y (binding feature). In this mode and bind.assay is set, normlised counts from layer 'data' will be used for comparing.  If bind.assay is none, will use raw counts from layer 'counts'. Mode 2: Compares X vs (Y - X). Mode 3: Compares X vs (Y + X). For mode 2 and 3, only support to use raw counts from layer 'counts'.
#' @param library.size Library size for each cell, used for normalise counts. If not set, use colSum(counts) instead.
#' @param wm.name Weight matrix name, this matrix (graph) generated by \code{\link{RunAutoCorr}}.
#' @param perm Permutations for evaluating mean and sd of D/L scores. Default is 100.
#' @param seed Seed for generate random number. Default is 999.
#' @param threads Threads. If set to 0 (default), will auto check the CPU cores and set threads = number of CPU cores -1.
#' @param versbose Print log message. Default is TRUE.
#' @param debug Print debug message. Will auto set thread to 1. Default is FALSE.
#' @param force Force to use this function.
#' @importFrom Matrix sparseMatrix
#' @export
RunBlockCorr <- function(object = NULL,
                         bind.name = "gene_name",
                         features = NULL,
                         assay = NULL,
                         wm.name = NULL,
                         bind.assay = NULL,
                         bind.features = NULL,
                         min.cells.bind = 10,
                         prefix = NULL,
                         min.features.per.block = 1,
                         scale.factor = 1e4,
                         mode = c(1,2,3),
                         library.size = NULL,
                         perm=100,
                         seed=999,
                         threads = 0,
                         verbose = TRUE,
                         debug = FALSE,
                         force = FALSE
                         )
{
  if (!force) {
    stop("RunBlockCorr is deprecated, please use RunSDT instead.")
  }
  mode <- mode[1L]
  tt <- Sys.time()
  
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  if (isTRUE(verbose)) {
    message("Working on assay ", assay, ".")
  }

  if (!is.null(bind.assay)) {
    if (bind.assay %ni% names(object)) {
      stop("No bind.assay is found. Make sure you specify the correct assay name.")
    } else {
      if (isTRUE(verbose)) {
        message("Working on binding assay ", bind.assay, ".")
      }
    }
  }
  
  wm.name <- wm.name %||% grep("_wm$", names(object), value = TRUE)
  if (length(wm.name) == 0) {
    stop("No weight matrix found. Perform RunAutoCorr() first.")
  }
  wm.name <- wm.name[1L]
  if (wm.name %ni% names(object)) {
    stop("No weight matrix found. Perform RunAutoCorr() first.")
  }
  if (verbose) {
    message("Use predefined weight matrix \"", wm.name, "\"", ".")
  }
  W <- object[[wm.name]]
  
  ## cell order may be changed due to merge, here reorder them
  cells <- colnames(W)
  
  features <- features %||% AutoCorrFeatures(object)
  features <- intersect(features, rownames(object))

  if (isTRUE(verbose)) {
    message("Processing ", length(features), " features.")
  }
  threads <- getCores(threads)
  
  object0 <- object[[assay]]
  tab <- object0[[]]

  if (bind.name %ni% colnames(tab)) {
    stop("No bind.name found in the feature table of assay ", assay, ". Run ParseExonName or ParseVarName first.")
  }

  bind.assay <- bind.assay %||% "tmp.assay"
  
  if (bind.assay %ni% names(object)) {
    if (min.features.per.block == 1) {
      if (isTRUE(verbose)) {
        message("No bind.assay specified, update min.features.per.block to 2.")
      }
      min.features.per.block <- 2
    }
  }
  # skip unannotated records
  tab <- tab[tab[[bind.name]] != "." & !is.na(tab[[bind.name]]),] 
  blocks <- names(which(table(tab[[bind.name]]) >= min.features.per.block))
  
  bind.features <- bind.features %||% blocks
  bind.features <- intersect(bind.features, blocks)
  bind.features <- intersect(bind.features, tab[features, bind.name])
  tab <- tab[which(tab[[bind.name]] %in% bind.features),]
  features <- intersect(features, rownames(tab))

  if (length(features) == 0) {
    stop("No features found.")
  } 
  
  if (isTRUE(verbose)) {
    message("Processing ", length(bind.features), " binding features.")
  }

  cs <- library.size %||% colSums(object, slot = "counts")

  norm <- TRUE
  if (bind.assay %ni% names(object)) {
    if (isTRUE(verbose)) {
      message("Use \"counts\" layer for test features.")
      message("Aggregate counts for binding features.")
    }
    x <- GetAssayData1(object, assay = assay, layer = "counts")
    features0 <- rownames(tab)
    x0 <- x[features0, cells]
    x0 <- as(x0, "TsparseMatrix")
    
    # Aggregate features in the same block
    y <- sparseMatrix(i = match(tab[[bind.name]][x0@i+1], bind.features),
                      j = x0@j+1,
                      x = x0@x, dims=c(length(bind.features), length(cells)))

    rownames(y) <- bind.features
    colnames(y) <- cells
    rm(x0)
    x <- x[features, cells]
  } else {
    if (isTRUE(verbose)) {
      message("Retrieve binding data from assay ", bind.assay, ".")
    }
    old.assay <- DefaultAssay(object)
    DefaultAssay(object) <- bind.assay

    ## todo
    y <- GetAssayData1(object, assay = bind.assay, layer = "counts")
    y <- y[,cells]
    rs <- Matrix::rowSums(y>0)
    idx <- which(rs >= min.cells.bind)
    blocks1 <- rownames(object)[idx]

    tab <- subset(tab, tab[[bind.name]] %in% blocks1)
    
    DefaultAssay(object) <- old.assay

    features <- intersect(features, rownames(tab))
    tab <- tab[features, ]
    bind.features <- unique(tab[[bind.name]])
    
    if (mode == 1) {
      if (isTRUE(verbose)) {
        message("Use \"data\" layer for test features and binding features.")
      }
      x <- GetAssayData1(object, assay = assay, layer = "data")
      y <- GetAssayData1(object, assay = bind.assay, layer = "data")
      x <- x[features, cells]
      y <- y[bind.features, cells]
      norm <- FALSE
    } else {
      x <- GetAssayData1(object, assay = assay, layer = "counts")  
      x <- x[features, cells]
    }
  }

  tab <- tab[features,]
  
  bidx <- match(tab[[bind.name]],rownames(y))
  idx <- match(features, rownames(x))
  if (isTRUE(verbose)) {
    message("Using ", threads, " threads.")
  }
  gc()

  tab <- object0[[]]
  
  ta <- .Call("D_test_v1", x, y, W, 1, perm, threads, idx, bidx, cs[cells], scale.factor, mode, FALSE, norm, seed, debug)
  if (length(ta) == 1) stop(ta[[1]])
  
  e <- ta[[1]]
  tval <- ta[[2]]
  
  names(tval) <- features
  
  pval <- pt(tval, df = perm - 1, lower.tail = FALSE)
  names(pval) <- features
  prefix <- prefix %||% bind.name
  tab[[paste0(prefix, ".t")]] <- tval[rownames(object)]
  tab[[paste0(prefix, ".pval")]] <- pval[rownames(object)]
  tab[[paste0(prefix, ".padj")]] <- p.adjust(pval[rownames(object)], method = "BH")
  rm(ta)
  
  gc()
  object0[[colnames(tab)]] <- tab
  object[[assay]] <- object0
  features <- head(features)
  object <- LogSeuratCommand(object)

  DefaultAssay(object) <- old.assay
  tt <- Sys.time()-tt
  if (isTRUE(verbose)) {
    message("Runtime : ",format(tt), ".");
  }
  object
}
makeModeData <- function(y, x, idx, bidx, mode)
{
  y <- y[bidx,]
  if (mode == 2) {
    y <- y - x[idx,]
  } else if (mode == 3) {
    y <- y + x[idx,]
  }
  rownames(y) <- rownames(x)[idx]
  y
}

lognorm <- function(mat = NULL, cs = NULL, scale.factor = 1e5)
{
  if (is.null(mat)) {
    stop("No matrix.")
  }
  if (is.null(cs)) {
    cs <- colSums(mat)
  }

  if ("dgCMatrix" %ni% class(mat)) {
    mat <- as(mat, 'CsparseMatrix')
  }
  nor <- .Call("lognorm", mat, cs, scale.factor)
  colnames(nor) <- colnames(mat)
  rownames(nor) <- rownames(mat)
  nor
}

#' @title RunSDT
#' @description Run spatial dissimilarity test for features and their binding features in parallel.
#' @param object Seurat object
#' @param bind.name The title of binding features in the meta table is important, as most users begin Yano to perform alternative splicing analysis. By default, the `bind.name` is set to "gene_name".
#' @param features Vector of features to calculate. Default is AutoCorrFeatures(object).
#' @param assay Work assay.
#' @param bind.assay Name of binding assay. If the binding assay is not set, raw counts of features from the same block will be aggregated first, followed by normalization.
#' @param bind.features List of bind features. If set, bind features will be subset with this list.
#' @param min.cells.bind Binding features detected in few than minimal number of cells will be skipped. Default is 10.
#' @param prefix Prefix name for output scores and values. Default is same with bind.name.
## @param subset Rules for subsetting meta table before selecting features to perform test.
#' @param scale.factor The scale factor used to normalize counts is set to a default of 1e4. The raw counts from "counts" layer are further normalized by sample size and scale factor for the spatial dissimilarity test.
#' @param mode Test mode. Mode 1: Compares X (test feature) vs Y (binding feature). In this mode and bind.assay is set, normlised counts from layer 'data' will be used for comparing.  If bind.assay is none, will use raw counts from layer 'counts'. Mode 2: Compares X vs (Y - X). Mode 3: Compares X vs (Y + X). For mode 2 and 3, only support to use raw counts from layer 'counts'.
#' @param library.size Library size for each cell, used for normalise counts. If not set, use colSum(counts) instead.
#' @param wm.name Weight matrix name, this matrix (graph) generated by \code{\link{RunAutoCorr}}.
#' @param perm Permutations for evaluating mean and sd of D/L scores. Default is 100.
#' @param seed Seed for generate random number. Default is 999.
#' @param threads Threads. If set to 0 (default), will auto check the CPU cores and set threads = number of CPU cores -1.
#' @param versbose Print log message. Default is TRUE.
#' @param debug Print debug message. Will auto set thread to 1. Default is FALSE.
#' @examples
#' data("glbt_small")
#' DefaultAssay(glbt_small) <- "RNA"
#' glbt_small <- NormalizeData(glbt_small) %>% RunUMAP(dim = 1:20)
#' DefaultAssay(glbt_small) <- "exon"
#' glbt_small <- NormalizeData(glbt_small)
#' glbt_small <- ParseExonName(glbt_small)
#' glbt_small <- RunAutoCorr(glbt_small)
#' glbt_small <- RunSDT(glbt_small, bind.name = "gene_name", bind.assay = "RNA")
#'
#' @importFrom Matrix sparseMatrix
#' @export
RunSDT <- function(object = NULL,
                   bind.name = "gene_name",
                   features = NULL,
                   assay = NULL,
                   wm.name = NULL,
                   bind.assay = NULL,
                   bind.features = NULL,
                   min.cells.bind = 10,
                   prefix = NULL,
                   min.features.per.block = 1,
                   scale.factor = 1e4,
                   mode = c(1,2,3),
                   library.size = NULL,
                   perm=100,
                   seed=999,
                   threads = 0,
                   verbose = TRUE,
                   debug = FALSE
                   )
{
  mode <- mode[1L]
  tt <- Sys.time()
  
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  if (isTRUE(verbose)) {
    message("Working on assay ", assay, ".")
  }

  if (!is.null(bind.assay)) {
    if (bind.assay %ni% names(object)) {
      stop("No bind.assay is found. Make sure you specify the correct assay name.")
    } else {
      if (isTRUE(verbose)) {
        message("Working on binding assay ", bind.assay, ".")
      }
    }
  }
  
  wm.name <- wm.name %||% grep("_wm$", names(object), value = TRUE)
  if (length(wm.name) == 0) {
    stop("No weight matrix found. Perform RunAutoCorr() first.")
  }
  wm.name <- wm.name[1L]
  if (wm.name %ni% names(object)) {
    stop("No weight matrix found. Perform RunAutoCorr() first.")
  }
  if (verbose) {
    message("Use predefined weight matrix \"", wm.name, "\"", ".")
  }
  W <- object[[wm.name]]
  
  ## cell order may be changed due to merge, here reorder them
  cells <- colnames(W)
  
  features <- features %||% AutoCorrFeatures(object)
  features <- intersect(features, rownames(object))

  if (isTRUE(verbose)) {
    message("Processing ", length(features), " features.")
  }
  threads <- getCores(threads)
  
  object0 <- object[[assay]]
  tab <- object0[[]]

  if (bind.name %ni% colnames(tab)) {
    stop("No bind.name found in the feature table of assay ", assay, ". Run ParseExonName or ParseVarName first.")
  }

  bind.assay <- bind.assay %||% "tmp.assay"
  
  if (bind.assay %ni% names(object)) {
    if (min.features.per.block == 1) {
      if (isTRUE(verbose)) {
        message("No bind.assay specified, update min.features.per.block to 2.")
      }
      min.features.per.block <- 2
    }
  }
  # skip unannotated records
  tab <- tab[tab[[bind.name]] != "." & !is.na(tab[[bind.name]]),] 
  blocks <- names(which(table(tab[[bind.name]]) >= min.features.per.block))
  blocks <- intersect(blocks, tab[features, bind.name])
  bind.features <- bind.features %||% blocks
  bind.features <- intersect(bind.features, blocks)
  tab <- tab[which(tab[[bind.name]] %in% bind.features),]
  features <- intersect(features, rownames(tab))

  if (length(features) == 0) {
    stop("No features found.")
  }
  
  if (isTRUE(verbose)) {
    message("Processing ", length(bind.features), " binding features.")
  }

  if (bind.assay %ni% names(object)) {
    if (isTRUE(verbose)) {
      message("Use \"counts\" layer for test features.")
      message("Aggregate counts for binding features.")
    }
    x <- GetAssayData1(object, assay = assay, layer = "counts")
    features0 <- rownames(tab)
    x0 <- x[features0, cells]
    x0 <- as(x0, "TsparseMatrix")
    
    # Aggregate features in the same block
    y <- sparseMatrix(i = match(tab[[bind.name]][x0@i+1], bind.features),
                      j = x0@j+1,
                      x = x0@x, dims=c(length(bind.features), length(cells)))

    rownames(y) <- bind.features
    colnames(y) <- cells
    rm(x0)
    if (mode == 1) {
      x <- GetAssayData1(object, assay = assay, layer = "data")
      x <- x[features, cells]
      cs <- library.size %||% colSums(object, slot = "counts")
      y0 <- lognorm(y, cs[cells], scale.factor)
      rm(y)
      y <- y0
    }
  } else {
    if (isTRUE(verbose)) {
      message("Retrieve binding data from assay ", bind.assay, ".")
    }
    old.assay <- DefaultAssay(object)
    DefaultAssay(object) <- bind.assay

    y <- GetAssayData1(object, assay = bind.assay, layer = "counts")
    y <- y[,cells]
    rs <- Matrix::rowSums(y>0)
    idx <- which(rs >= min.cells.bind)
    blocks1 <- rownames(y)[idx]

    tab <- subset(tab, tab[[bind.name]] %in% blocks1)
    
    DefaultAssay(object) <- old.assay

    features <- intersect(features, rownames(tab))
    tab <- tab[features,]
    bind.features <- unique(tab[[bind.name]])
    
    if (mode == 1) {
      if (isTRUE(verbose)) {
        message("Use \"data\" layer for test features and binding features.")
      }
      x <- GetAssayData1(object, assay = assay, layer = "data")
      y <- GetAssayData1(object, assay = bind.assay, layer = "data")
      x <- x[features, cells]
      y <- y[bind.features, cells]
    } else {
      x <- GetAssayData1(object, assay = assay, layer = "counts")
      x <- x[features, cells]
    }
  }

  ## only keep features in the table
  tab <- tab[features,]
  
  bidx <- match(tab[[bind.name]],rownames(y))
  idx <- match(features, rownames(x))
  if (mode != 1) {
    cs <- library.size %||% colSums(object, slot = "counts")      
    y <- makeModeData(y, x, idx, bidx, mode)
    y0 <- lognorm(y, cs[cells], scale.factor)
    rm(y)
    y <- y0
    x <- lognorm(x[idx,], cs[cells], scale.factor)
    idx <- 1:nrow(x)
    bidx <- 1:nrow(y)
  }
  
  if (isTRUE(verbose)) {
    message("Using ", threads, " threads.")
  }
  gc()
  
  tab <- object0[[]]
  
  ta <- .Call("D_test_v2", x, y, W, perm, threads, idx, bidx, 0, seed, debug)
  if (length(ta) == 1) stop(ta[[1]])
  
  d <- ta[[1]]
  tval <- ta[[2]]
  names(d) <- features
  names(tval) <- features
  
  pval <- pt(tval, df = perm - 1, lower.tail = FALSE)
  names(pval) <- features
  prefix <- prefix %||% bind.name
  tab[[paste0(prefix, ".D")]] <- d[rownames(object)]
  tab[[paste0(prefix, ".t")]] <- tval[rownames(object)]
  tab[[paste0(prefix, ".pval")]] <- pval[rownames(object)]
  tab[[paste0(prefix, ".padj")]] <- p.adjust(pval[rownames(object)], method = "BH")
  rm(ta)
  
  gc()
  object0[[colnames(tab)]] <- tab
  object[[assay]] <- object0
  features <- head(features)
  object <- LogSeuratCommand(object)

  DefaultAssay(object) <- old.assay
  tt <- Sys.time()-tt
  if (isTRUE(verbose)) {
    message("Runtime : ",format(tt), ".");
  }
  object
}

cor_dist <- function(x = NULL, y = NULL, W = NULL, perm = 1000, thread = 1)
{
  ta <- .Call("D_distribution_test", x, y, W, perm, thread)
  ta
}
cor_dist2 <- function(x = NULL, y = NULL, W = NULL, perm = 1000, thread = 1)
{
  ta <- .Call("D_distribution_test2", x, y, W, perm, thread)
  ta
}

#' @export
SDTdemo <- function(object = NULL, bind.name = NULL, bind.assay = NULL, assay = NULL, mode = c(1,2,3), perm = 100, feature = NULL, library.size = NULL, method = c("D1","D2")) {
  
  if (is.null(object)) {
    stop("No object.")
  }

  if (is.null(feature)) {
    stop("No feature.")
  }

  if (is.null(bind.name)) {
    stop("bind.name is not set.")
  }
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)

  mode <- mode[1L]
  method <- method[1L]
  
  DefaultAssay(object) <- assay
  feature <- intersect(rownames(object), feature)
  feature <- feature[1L]
  
  if (length(feature) == 0) {
    stop("No feature found.")
  }

  df <- object[[assay]][[]]

  if (!(bind.name %in% colnames(df))) {
    stop("No bind.name found.")
  }

  bind.feature <- df[feature, bind.name]

  if (is.na(bind.feature)) {
    stop("bind.feature is NA.")
  }

  if (mode == 1) {
    x <- GetAssayData(object, layer = "data")
    DefaultAssay(object) <- bind.assay
    y <- GetAssayData(object, layer = "data")
  } else {
    x <- GetAssayData(object, layer = "counts")
    DefaultAssay(object) <- bind.assay
    y <- GetAssayData(object, layer = "counts")
    cs <- library.size %||% colSums(x)
  }
  
  d1 <- x[feature,]
  d1 <- as.matrix(d1)

  bind.feature <- intersect(bind.feature, rownames(object))
  if (length(bind.feature) == 0) {
    stop("No bind.feature found.")
  }

  d2 <- y[bind.feature,]
  d2 <- as.matrix(d2)
  if (mode ==  2) {
    d2 <- d2 - d1
  }
  if (mode == 3) {
    d2 <- d1 + d2
  }
  
  message("Orginal cor is ", cor(d1[,1], d2[,1]))

  if (mode != 1) {
    d1[,1] <- log(d1[,1]*1e4/cs+1)
    d2[,1] <- log(d2[,1]*1e4/cs+1)
  }
  
  message("Cor after normlise is ", cor(d1[,1], d2[,1]))
  message("Mean a ", mean(d1[,1]), " mean b ", mean(d2[,1]))
  
  W <- object[['pca_wm']]
  W <- as(W, "CsparseMatrix")

  s1 <- (d1[,1] %*% W)[1,]
  s2 <- (d2[,1] %*% W)[1,]
  
  message("Cor after smooth is ", cor(s1, s2))
  message("Mean smooth a ", mean(s1), " mean smooth b ", mean(s2))
  d1 <- d1[,1]
  d2 <- d2[,1]
  sm <- mean(d1)
  sm2 <- mean(d2)
  Lx <- sum((s1-sm)^2)/sum((d1-sm)^2)
  Ly <- sum((s2-sm2)^2)/sum((d2-sm2)^2)

  D <- sqrt(Lx)*(1-cor(s1,s2))
  D2 <- sqrt(Ly)*(1-cor(s1,s2))
  
  mn <- lapply(1:perm, function(i) {
    d11 <- sample(d1)
    s11 <- (d11 %*% W)[1,]
    Lx1 <- sum((s11-sm)^2)/sum((d11-sm)^2)
    sqrt(Lx1)*(1-cor(s11,s2))    
  })
  
  mn <- unlist(mn)
  hist(mn)
  m <- mean(mn)
  var <- sqrt(sum((mn - m)^2)/perm)
  t <- (D-m)/var
  
  p <- pt(t, df = perm-1, lower.tail = FALSE)
  message("Lx is ", Lx, ", D score is ", D)
  message("Mean is ", m, ", var is ", var)
  message("t value is ", t, ";\np value is ", p)

  if (method == "D2") {
    message("================== D2 ==========")
    mn <- lapply(1:perm, function(i) {
      d22 <- sample(d2)
      s22 <- (d22 %*% W)[1,]
      Ly1 <- sum((s22-sm2)^2)/sum((d22-sm2)^2)
      sqrt(Ly1)*(1-cor(s1,s22))    
    })

    mn <- unlist(mn)
    hist(mn)
    m <- mean(mn)
    var <- sqrt(sum((mn - m)^2)/perm)
    t <- (D2-m)/var
  
    p <- pt(t, df = perm-1, lower.tail = FALSE)
    message("Ly is ", Ly, ", D2 score is ", D2)
    message("Mean is ", m, ", var is ", var)
    message("t value is ", t, ";\np value of D2 is ", p)
  }
  DefaultAssay(object) <- old.assay
}

