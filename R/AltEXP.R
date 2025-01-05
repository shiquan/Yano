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
  min.cells = 0
) {
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) +length(x=cells.2)< min.cells) {
    stop("Cell group 1 has fewer than ", min.cells, " cells")
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

DEXSeqTest <- function(x, y, cells.1 = NULL, cells.2 = NULL, pseudo.group = 3, rst = NULL, mode = 1, threads = 1, ...)
{
  if (!PackageCheck('DEXSeq', error = FALSE)) {
    stop("Install DEXSeq package first. BiocManager::install('DEXSeq')")
  }
  
  stopifnot(identical(dim(x), dim(y)))

  n1 <- length(cells.1)
  n2 <- length(cells.2)

  new.group1 <- rep(1:pseudo.group,(n1/pseudo.group+1), length.out = n1)
  new.group2 <- rep((pseudo.group+1):(pseudo.group*2),(n2/pseudo.group+1), length.out = n2)
  
  new.group <- c(new.group1,new.group2)

  features <- rownames(x)

  x <- x[, c(cells.1, cells.2)]
  x <- as(x,"TsparseMatrix")
  x <- Matrix::sparseMatrix(i=(x@i+1),j=new.group[x@j+1], x=x@x)
  colnames(x) <- paste0("group_",1:(pseudo.group*2))

  y <- y[, c(cells.1, cells.2)]
  y <- as(y, "TsparseMatrix")
  y <- Matrix::sparseMatrix(i=(y@i+1),j=new.group[y@j+1], x=y@x)
  colnames(y) <- paste0("group_",1:(pseudo.group*2))

  meta.tab <- data.frame(row.names=colnames(x),
                         condition=factor(c(rep("cluster1", pseudo.group),
                                            rep("cluster2", pseudo.group))))
  
  if (mode == 2) {
    y <- y - x
    y[y<0] <- 0
  }

  if (mode == 3) {
    y <- y + x
  }

  meta.tab$condition <- factor(meta.tab$condition)

  feature.names <- paste0("feature_",1:nrow(x))
  names(feature.names) <- features
  rownames(x) <- feature.names
  rownames(y) <- feature.names

  dxd <- DEXSeq::DEXSeqDataSet(as.matrix(x),
                               sampleData=meta.tab,
                               groupID = rst$bind.feature,
                               featureID = feature.names,
                               alternativeCountData = as.matrix(y),
                               design= ~sample+exon+condition:exon)

  rm(x)
  rm(y)

  # threads <- getCores(threads)
  # BPPARAM <- BiocParallel::MulticoreParam(workers = threads)

  results <- DEXSeq::DEXSeq(dxd, ...)
  if (is.character(results)) {
    return(results)
  }
  rst$log2fc <- results$log2fold_cluster2_cluster1
  rst$pval <- results$pvalue
  rst$padj <- results$padj

  return(rst)
}
PermTest <- function(x, y, cells.1, cells.2, rst, perm = 100, seed = 999, mode = c(1,2,3), threads = 1, debug = FALSE)
{
  stopifnot(identical(dim(x), dim(y)))
  x <- as(x, "dgCMatrix")
  y <- as(y, "dgCMatrix")

  cells <- colnames(x)
  if (is.null(cells)) {
   stop("Colnames of X is empty.")
  }
  
  idx.1 <- match(cells.1, cells)
  idx.2 <- match(cells.2, cells)
  
  df <- .Call("alt_exp", x, y, idx.1, idx.2, mode[1L], perm, threads, seed)

  if (length(df) == 1) {
    stop(df)
  }
  rst$delta <- df[[1]]
  rst$tval <- df[[2]]
  rst$mean <- df[[3]]
  rst$var <- df[[4]]
  rst$ratio.1 <- df[[5]]
  rst$ratio.2 <- df[[6]]
  rst$pval <- 2*(1-pt(abs(rst$tval), df = perm - 1))
  rst$padj <- p.adjust(rst$pval, method = "BH")
  return (rst)
}

#'
#' @name RunDEXSeq
#' @title RunDEXSeq
#' @description This method tests whether the test features and binding features are differentially expressed between groups using a generalized linear model. For details, refer to the DEXSeq method (PMID: 22722343). Since this approach is applied to single cells, cells within the same group are aggregated into several pseudo-bulk datasets for analysis.
#' @param object A Seurat object.
#' @param bind.name The title of binding name in meta table. Usually be "gene_name" for alternative splicing.
#' @param ident.1 Identify class to test, if not set will compare all groups one by one
#' @param ident.2 A second class for comparsion. If NULL (default), use all other cells for comparison.
#' @param cells.1 Vector of cell names belong to group 1. Conflict with ident.1
#' @param cells.2 Vector of cell names for comparsion. Conflict with ident.2
#' @param assay Assay for test features. Default assay will be used if not set.
#' @param bind.assay Assay for binding features. If not set, test features in same goup (with same bind name) will be aggreated as binding feature
#' @param features Candidate list to test. If not set, will use all.
#' @param bind.features Candidate list for binding features. If not set, will use all.
#' @param min.pct Only test features that are detected in a minimum fraction of min.pct cells in all cells.  Meant to speed up the function by not testing genes that are very infrequenctly expressed in all cells. Remember we are testing alternative epxression pattern here, so it is possible the test feature is not expressed in one group, therefore we are not going to check by groups. Note that min.pct is set for test feature here. But in \code{\link{RunPSI}}, the min.pct is set for binding feature. Default is 0.05.
#' @param min.pct.bind.feature Only test binding features that are detected in a minimum fraction of min.pct.bind.feature in either of the two populations. Meant to speed up the function by not testing genes that are very infrequenctly expressed in both groups. Default is 0.05.
#' @param return.thresh Only return markers that have a p-value < return.thresh.
#' @param node A node to find markers for and all its children; requires \code{\link{BuildClusterTree}} to have been run previously. Only can be used if test all groups.
#' @param pseudo.group Aggregate single cells into pseudo groups, because DEXSeq is designed for bulk RNA-seq. At least 3 cells are required for each group. Default is 3.
#' @param mode Test mode, default is 1. See online manual for the difference between modes. <https://shiquan.github.io/Yano.html>
#' @param threads Threads passed to DEXSeq. Default is 1.
#' @return Data frame containing p values and pct for test features and their binding features.
#' @export
RunDEXSeq <- function(object = NULL, bind.name = "gene_name",
                      ident.1 = NULL, ident.2 = NULL, cells.1 = NULL, cells.2 = NULL,
                      assay = NULL, bind.assay = NULL,
                      features = NULL, bind.features = NULL,
                      min.pct = 0.05, min.pct.bind.feature = 0.05,
                      return.thresh = NULL, node = NULL,
                      pseudo.group = 3, mode = 1, threads = 1, debug = FALSE)
{
  if (is.null(object)) {
    stop("No object specified.")
  }

  if (is.null(bind.name)) {
    stop("No bind.name specified.")
  }

  if (is.null(bind.assay)) {
    warnings("Bind assay is not specified, will aggreate exons from the same gene.")
  }
  
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  
  df <- object[[assay]][[]]
  if (bind.name %ni% colnames(df)) {
    stop("No bind.name found in meta data. You may need ParseExonName first.")
  }

  if (!is.null(ident.1) | !is.null(cells.1)) {
    tb <- FindAltExp(object, ident.1 = ident.1, ident.2 = ident.2, cells.1 = cells.1, cells.2 = cells.2,
                     assay = assay, bind.assay = bind.assay,
                     bind.name = bind.name,
                     test.use = "DEXSeq",
                     #min.pct = min.pct, min.pct.bind.feature = min.pct.bind.feature,
                     mode = mode,
                     features = features,
                     return.thresh = return.thresh,
                     pseudo.group = pseudo.group, threads = threads,
                     debug = debug)
  } else {
    tb <- FindAllAltExp(object, assay = assay, bind.assay = bind.assay, bind.name = bind.name,
                        test.use = "DEXSeq",
                        node = node, features = features,
                        return.thresh = return.thresh,
                        #min.pct = min.pct, min.pct.bind.feature = min.pct.bind.feature,
                        mode = mode,
                        pseudo.group = pseudo.group, threads = threads, debug = debug)
  }
  tb
}

#' @name FindAltExp
#' @title Find Alternative expression
#' @description Using linear regression model (implemented by DEXSeq) or PermTest method to calculate delta ratio and simulate p value for each event.
#' @details This function first aggregates all raw counts per feature for each group and then perform alternative expression analysis with DEXSeq or PermTest method. For DEXSeq, counts per group will divide into pesudo-groups first. For PermTest method, \eqn{delta ratio = X_a/Y_a - X_b/Y_b}. It then employs a permutation method to randomize the cells in the two groups 100 times (by default) to evaluate the mean and standard deviation of delta-ratio. A p-value is calculated using a t-test. The PermTest method is not well tested. Therefore, we use DEXSeq method in default, but it's very slow if input a lot of features.
#' @param object A Seurat object.
#' @param ident.1 Identify class to test, if not set will compare all groups one by one
#' @param ident.2 A second class for comparsion. If NULL (default), use all other cells for comparison.
#' @param cells.1 Vector of cell names belong to group 1. Conflict with ident.1
#' @param cells.2 Vector of cell names for comparsion. Conflict with ident.2
#' @param assay Test assay (X). Default assay will be used if not set.
#' @param bind.name Title name for binding features in the meta table. Consider most users start Yano to perform alternative splicing analysis, the default bind.name set to "gene_name".
#' @param bind.assay Bind assay (Y). If not set, will aggregate all X values of the same block.
#' @param features Candidate list to test. If not set, will use AutoCorrFeatures(object, assay = assay).
#' @param bind.features Candidate list for bind features to test. If not set, will test all covered.
#' @param min.cells Used to filter candiate features or binding features. Require them at least expressed in min.cells. Default is 10.
#' @param pesudo.group Aggregate counts into groups for each clusters. Used only for DEXSeq.
#' @param return.thresh Only return markers that have a p-value < return.thresh. Default is NULL.
#' @param mode Test mode. For mode 1, X (test feature) vs Y (binding feature). For mode 2, X vs (Y-X). For mode 3, X vs (Y+X). 
#' @param perm Permutation steps for calculate statistical of delta-ratio. Default is 100.
#' @param seed Seed for generate random number. Default is 999.
#' @param threads Threads. For DEXSeq, threads will set to 1. For other methods, threads set to 0, which will auto check the CPU cores and set threads = number of CPU cores -1.
#' @param debug Print debug logs. Will auto set thread to 1. Default is FALSE.
#' @importFrom SeuratObject PackageCheck
#' @export
FindAltExp <- function(object = NULL,
                       cells.1 = NULL, cells.2 = NULL,
                       ident.1 = NULL, ident.2 = NULL,
                       assay = NULL,
                       bind.name = "gene_name",
                       test.use = c("DEXSeq", "PermTest"),
                       bind.assay = NULL,
                       features = NULL,
                       bind.features = NULL,
                       min.cells = 10,
                       return.thresh = NULL,
                       mode = c(1,2,3),
                       threads = 0,
                       perm = 100,
                       seed = 999,
                       pseudo.group = 3,
                       debug = FALSE
                       )
{
  if (is.null(bind.name)) {
    stop("No bind.name specified.")
  }

  test.use <- match.arg(test.use)
  
  assay <- assay %||% DefaultAssay(object)
  
  bind.features0 <- CheckBindName(object, bind.name, assay = assay)
  
  if (!is.null(ident.1) & !is.null(cells.1)) {
    stop("cells.1 is conflict with ident.1.")
  }

  if (!is.null(ident.2) & !is.null(cells.2)) {
    stop("cells.2 is conflict with ident.2.")
  }

  if (!is.null(ident.1)) {
    cells.1 <- colnames(object)[which(Idents(object) == ident.1)]
    if (is.null(cells.1)) stop("No ident.1 found")
  }

  if (!is.null(ident.2)) {
    cells.2 <- colnames(object)[which(Idents(object) == ident.2)]
    if (is.null(cells.2)) stop("No ident.2 found")
  }

  if (is.null(cells.1)) {
    stop("No cells.1 or ident.1 defined.")
  }

  if (is.null(cells.2)) {
    cells.2 <- setdiff(colnames(object), cells.1)
  }
  
  ValidateCellGroups(object, cells.1, cells.2, min.cells)

  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  object0 <- object[[assay]]
  df <- object0[[]]
  if (bind.name %ni% colnames(df)) {
    stop("No bind name found in meta data. Run ParseExonName or ParseVarName first.")
  }
  
  features <- features %||% AutoCorrFeatures(object, assay = assay)
  features <- intersect(features, rownames(object0))
  # filter features
  dat <- GetAssayData1(object, layer="counts")
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
  rst$bind.feature <- bind.features0[features]
  # rst <- rst[which(pct >= min.pct & !is.na(rst$bind.feature)),]
  
  if (!is.null(bind.features)) {
    rst <- subset(rst, bind.feature %in% bind.features)
  }
  
  features <- rownames(rst)
  bind.features <- unique(rst$bind.feature)
  
  x <- x[features, ]
  
  if (is.null(bind.assay)) {

    idx <- which(df[[bind.name]] %in% bind.features)
    df0 <- df[idx,]
    idx <- names(which(table(df0[[bind.name]])>1))
    idx <- which(df[[bind.name]] %in% idx)
    df0 <- df[idx,]
    features0 <- rownames(df0)
    x0 <- dat[features0,]
    x0 <- as(x0, "TsparseMatrix")
    arr <- df0[[bind.name]]
    bind.features <- unique(arr)
    # Aggregate features in the same block
    y <- sparseMatrix(i = match(arr[x0@i+1], bind.features),
                      j = x0@j+1,
                      x = x0@x, dims=c(length(bind.features), ncol(dat)))

    rm(x0)
    colnames(y) <- colnames(dat)
    rownames(y) <- bind.features
    y <- y[, c(cells.1, cells.2)]
  } else {
    y <- GetAssayData1(object,assay = bind.assay, layer = "counts")
    y <- y[, c(cells.1, cells.2)]
    bind.features <- intersect(unique(bind.features), rownames(y))
  }
  rst <- subset(rst, bind.feature %in% bind.features)
  features <- rownames(rst)
  x <- x[features,]
  y <- y[rst$bind.feature,]
  stopifnot(identical(dim(x), dim(y)))

  mode <- mode[1L]
  
  if (mode == 2) {
    y <- y - x
    y[y<0] <- 0
  } else if (mode == 3) {
    y <- y + x
  }

  rownames(y) <- rownames(x)
  
  y <- y[rst$feature,]
  x <- x[rst$feature,]

  DefaultAssay(object) <- old.assay
  threads <- getCores(threads)

  if (test.use == "PermTest") {
    rst <- PermTest(x, y, cells.1, cells.2, rst, perm = perm, seed = seed, mode = 1, threads = threads)
  } else {
    rst <- DEXSeqTest(x, y, cells.1, cells.2, rst = rst, mode = 1, pseudo.group = pseudo.group)
  }
  
  if (!is.null(return.thresh)) {
    rst <- subset(rst, padj < return.thresh)
  }
  rst
}

#' @name FindAllAltExp
#' @title Test alternative expression for all cell groups
#' @inheritParams FindAltExp
#' @param node A node to find markers for and all its children; requires \code{\link{BuildClusterTree}} to have been run previously. Only can be used if test all groups.
#' @return Data frame containing p values.
#' @export
#'
#' @examples
#' data("glbt_small")
#' DefaultAssay(glbt_small) <- "exon"
#' alt.exon <- FindAllAltExp(object = glbt_small, bind.assay = "RNA", bind.name = "gene_name", features = rownames(glbt_small))
#' head(alt.exon)
#' 
#' @importFrom SeuratObject PackageCheck
#' @importFrom Seurat Tool
#' @export
FindAllAltExp <- function(object = NULL,
                          assay = NULL,
                          features = NULL,
                          bind.name = "gene_name",
                          bind.assay = NULL,
                          bind.features = NULL,
                          node = NULL,
                          min.cells = 10,
                          return.thresh = NULL,
                          mode = c(1,2,3),
                          test.use = c("DEXSeq","PermTest"),
                          threads = 0,
                          perm = 100,
                          seed = 999,
                          pseudo.group = 3,
                          debug = FALSE
                          )
{
  ## node paramenter inhert from Seurat v5.0.0, check Seurat::FindAllMarkers for details
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
    if (!PackageCheck('ape', error = FALSE)) {
      stop("Install ape package")
    }
    tree <- Seurat::Tool(object = object, slot = 'BuildClusterTree')
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
        FindAltExp(
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
          assay = assay,
          bind.name = bind.name,
          bind.assay = bind.assay,
          features = features,
          bind.features = bind.features,
          min.cells = min.cells,
          return.thresh = return.thresh,
          mode = mode[1L],
          test.use = test.use,
          threads = threads,
          pseudo.group = pseudo.group,
          perm = perm,
          seed = seed,
          debug = debug
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

  if (length(x = messages) > 0) {
    warning("The following tests were not performed: ", call. = FALSE, immediate. = TRUE)
    for (i in 1:length(x = messages)) {
      if (!is.null(x = messages[[i]])) {
        warning("When testing ", idents.all[i], " : ", messages[[i]], call. = FALSE, immediate. = TRUE)
      }
    }
  }

  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
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
