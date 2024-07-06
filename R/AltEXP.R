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
  return (rst)
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

  y <- y[, c(cells.1, cells.2)]
  y <- as(y, "TsparseMatrix")
  y <- Matrix::sparseMatrix(i=(y@i+1),j=new.group[y@j+1], x=y@x)
  colnames(y) <- paste0("group_",1:(pesudo.group*2))

  meta.tab <- data.frame(row.names=colnames(x),
                         condition=c(rep("cluster1", pesudo.group),
                                     rep("cluster2", pesudo.group)))
  
  if (length(x = features) == 0) {
    warning("No feature pass min.pct threshold; returning NULL.")
    return(NULL)
  }

  stopifnot(identical(dim(x), dim(y)))

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
AlternativeExpressionTest <- function(object,
                                      cells.1 = NULL, cells.2 = NULL,
                                      ident.1 = NULL, ident.2 = NULL,
                                      bind.name = NULL, bind.assay = NULL,
                                      assay = NULL,
                                      layer = "counts",
                                      test.use = c("DEXSeq", "PSI"),
                                      features = NULL,
                                      bind.features = NULL,
                                      min.pct = 0.01,
                                      min.cells.group = 3,
                                      min.pct.bind.feature = 0.01,
                                      pesudo.group = 3,
                                      return.thresh = 1e-2,
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

  if (!is.null(ident.1) & !is.null(cells.1)) {
    stop("cells.1 is conflict with ident.1")
  }

  if (!is.null(ident.2) & !is.null(cells.2)) {
    stop("cells.2 is conflict with ident.2")
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
  rst$bind.feature <- bind.features0[features]
  rst <- rst[which(pct >= min.pct & !is.na(rst$bind.feature)),]
  
  if (!is.null(bind.features)) {
    rst <- subset(rst, bind.feature %in% bind.features)
  }
  
  features <- rownames(rst)
  bind.features <- unique(rst$bind.feature)
  
  x <- x[features, ]
  
  if (is.null(bind.assay)) {
    x0 <- dat[features,]
    x0 <- as(x0, "TsparseMatrix")

    # Aggregate features in the same block
    y <- sparseMatrix(i = match(bind.features[x0@i+1], bind.features),
                      j = x0@j+1,
                      x = x0@x, dims=c(length(bind.features), ncol(dat)))

    rm(x0)
    colnames(y) <- colnames(dat)
    rownames(y) <- bind.features
    y <- y[, c(cells.1, cells.2)]
  } else {
    y <- GetAssayData1(object,assay = bind.assay, layer = layer)
    bind.features <- intersect(unique(bind.features), rownames(y))
    rst <- subset(rst, bind.feature %in% bind.features)
  }
  
  features <- rownames(rst)
  x <- x[features,]
  y <- y[rst$bind.feature,]
  stopifnot(identical(dim(x), dim(y)))
  
  if (mode == 2) {
    y <- y - x
    y[y<0] <- 0
  }

  if (mode == 3) {
    y <- y + x
  }

  rownames(y) <- rownames(x)
  b1 <- y[, cells.1]
  b2 <- y[, cells.2]

  bind.pct1 <- rowSums(b1>0)/length(cells.1)
  bind.pct2 <- rowSums(b2>0)/length(cells.2)
  rm(b1)
  rm(b2)
    
  names(bind.pct1) <- rownames(x)
  names(bind.pct2) <- rownames(y)
    
  idx <- which(bind.pct1 >= min.pct.bind.feature & bind.pct2 >= min.pct.bind.feature)
  rst <- rst[idx,]
  
  y <- y[rst$feature, c(cells.1, cells.2)]
  x <- x[rst$feature,]

  rst$bind.feature.pct.1 <- bind.pct1[rst$feature]
  rst$bind.feature.pct.2 <- bind.pct2[rst$feature]

  ## No matter which mode you set before, now y has been adjusted, so use mode ==1 here.
  if (method == "DEXSeq") {
    rst <- DEXSeqTest(x, y, cells.1, cells.2, pesudo.group, rst, mode = 1, threads = threads, ...)
  } else if (method == "PSI") {
    rst <- PermTest(x, y, cells.1, cells.2, rst, perm = perm, seed = seed, mode = 1, threads = threads, ...)
  }
  
  if (!is.null(return.thresh)) {
    rst <- subset(rst, padj < return.thresh)
  }
  rst
}

#'@importFrom SeuratObject PackageCheck
FindAllAEMarkers <- function(object,
                             bind.name = NULL, bind.assay = NULL,
                             layer = "counts",
                             features = NULL,
                             node = NULL,
                             min.pct = 0.01,
                             min.cells.group = 3,
                             min.pct.bind.feature = 0.01,
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

  ## node paramenter inhert from Seurat v5.0.0, check Seurat::FindAllMarkers for details
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

#' @name RunDEXSeq
#' @title Test for alternative expression with DEXSeq
#' @description This will test the test features and binding features are different expressed between groups with a generalised linear model. See DEXSeq and DESeq2 paper for details.
#' @param object A Seurat object.
#' @param bind.name The title of binding name in meta table.
#' @param ident.1 Identify class to test, if not set will compare all groups one by one
#' @param ident.2 A second class for comparsion. If NULL (default), use all other cells for comparison.
#' @param cells.1 Vector of cell names belong to group 1. Conflict with ident.1
#' @param cells.2 Vector of cell names for comparsion. Conflict with ident.2
#' @param assay Assay for test features. Default assay will be used if not set.
#' @param bind.assay Assay for binding features. If not set, test features in same goup (with same bind name) will be aggreated as binding feature
#' @param features Candidate list to test. If not set, will use all.
#' @param bind.features Candidate list for binding features. If not set, will use all.
#' @param min.pct Only test features that are detected in a minimum fraction of min.pct cells in all cells.  Meant to speed up the function by not testing genes that are very infrequenctly expressed in all cells. Remember we are testing alternative epxression pattern here, so it is possible the test feature is not expressed in one group, therefore we are not going to check by groups. Note that min.pct is set for test feature here. But in \code{\link{RunPSI}}, the min.pct is set for binding feature. Default is 0.01.
#' @param min.pct.bind.feature Only test binding features that are detected in a minimum fraction of min.pct.bind.feature in either of the two populations. Meant to speed up the function by not testing genes that are very infrequenctly expressed in both groups. Default is 0.01.
#' @param return.thresh Only return markers that have a p-value < return.thresh.
#' @param node A node to find markers for and all its children; requires \code{\link{BuildClusterTree}} to have been run previously. Only can be used if test all groups.
#' @param pesudo.group Aggregate single cells into pesudo groups, because DEXSeq is designed for bulk RNA-seq. At least 3 cells are required for each group. Default is 3.
#' @param mode Test mode, default is 1. See online manual for the difference between modes. <https://shiquan.github.io/Yano.html>
#' @param threads Threads passed to DEXSeq. Default is 1.
#' @return Data frame containing p values and pct for test features and their binding features.
#' @export
#'
#' @examples
#' data("neuron_small")
#' alt.exon <- RunDEXSeq(object = neuron_small, assay = "flatten", bind.assay = "RNA", bind.name = "gene_name")
#' head(alt.exon)
#' 
RunDEXSeq <- function(object = NULL, bind.name = "bind_name", ident.1 = NULL, ident.2 = NULL, cells.1 = NULL, cells.2 = NULL, assay = NULL, bind.assay = NULL, features = NULL, bind.features = NULL, min.pct = 0.01, min.pct.bind.feature = 0.01, return.thresh = 1e-2, node = NULL, pesudo.group = 3, mode = 1, threads = 1)
{
  if (is.null(object)) {
    stop("No object specified.")
  }

  if (is.null(bind.name)) {
    stop("No bind.name specified.")
  }
  
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  
  df <- object[[assay]][[]]
  if (bind.name %ni% colnames(df)) {
    stop("No bind.name found in meta data. You may need ParseExonName first.")
  }

  if (!is.null(ident.1) | !is.null(cells.1)) {
    tb <- AlternativeExpressionTest(object, ident.1 = ident.1, ident.2 = ident.2, cells.1 = cells.1, cells.2 = cells.2, assay = assay, bind.assay = bind.assay,
                                    bind.name = bind.name, test.use = "DEXSeq", min.pct = min.pct, min.pct.bind.feature = min.pct.bind.feature, mode = mode,
                                    return.thresh = return.thresh,
                                    pesudo.group = pesudo.group, threads = threads)
  } else {
    tb <- FindAllAEMarkers(object, assay = assay, bind.assay = bind.assay, bind.name = bind.name, test.use = "DEXSeq", node = node, features = features,
                           return.thresh = return.thresh,
                           min.pct = min.pct, min.pct.bind.feature = min.pct.bind.feature, mode = mode, pesudo.group = pesudo.group, threads = threads)
  }
  
  tb
}

#' @name RunPSI
#' @title Test for alternative expression with delta-PSI
#' @description This will calculate the delta-PSI = PSI_a - PSI_b for group a and b. Then use the permutation method to randomize the cells in two groups 100 times to evalue the mean and sd of delta-PSI and calculate the p value with t test. PSI is short for Percent-Spliced In, calculated by PSI = EXON/(EXON+EXCL), EXON indicates reads overlapped this exon, EXCL indicates excluded reads, which are junction reads skip this exon.
#' @param object A Seurat object.
#' @param ident.1 Identify class to test, if not set will compare all groups one by one
#' @param ident.2 A second class for comparsion. If NULL (default), use all other cells for comparison.
#' @param cells.1 Vector of cell names belong to group 1. Conflict with ident.1
#' @param cells.2 Vector of cell names for comparsion. Conflict with ident.2
#' @param exon.assay Assay for exons. Default assay will be used if not set.
#' @param exclude.assay Assay for exon exclude reads. This assay is mandantory to calculate PSI. The count matrix can be calculate with `PISA anno -exon -psi`.
#' @param features Candidate list to test. If not set, will use all.
#' @param genes Candidate list for genes to test. If not set, will test all covered.
#' @param exon.name Exon name in the title of meta table for exon assay. Default is "exon_name"
#' @param min.pct Excluded + exon reads that are detected in a minimum fraction of min.pct in either of the two populations. Meant to speed up the function by not testing genes that are very infrequenctly expressed in both groups. Note that the min.pct is for binding feature. Default is 0.01.
#' @param min.pct.exon Only test features that are detected in a minimum fraction of min.pct.exon cells in all cells.  Meant to speed up the function by not testing genes that are very infrequenctly expressed in all cells. Remember we are testing alternative epxression pattern here, so it is possible the exon is not expressed in one group, therefore we are not going to check by groups. Default is 0.01.
#' @param return.thresh Only return markers that have a p-value < return.thresh.
#' @param node A node to find markers for and all its children; requires \code{\link{BuildClusterTree}} to have been run previously. Only can be used if test all groups.
#' @param perm Permutation steps for calculate statistical of delta-PSI.
#' @param threads Threads to run. Default is 1.
#' @return Data frame containing p values and pct for exons and their genes.
#' @export
#'
#' @examples
#' data("neuron_small")
#' alt.exon <- RunPSI(object = neuron_small, exon.assay = "exon", exclude.assay = "EXCL")"
#' head(alt.exon)
#' 
RunPSI <- function(object = NULL, ident.1 = NULL, ident.2 = NULL, cells.1 = NULL, cells.2 = NULL, exon.assay = NULL, exclude.assay = NULL, features = NULL, genes = NULL, gene.name = "gene_name", exon.name = "exon_name", min.pct = 0.01, min.pct.exon = 0.01, return.thresh = 1e-2, node = NULL, perm = 100, seed = 999, threads = 1)
{
  if (is.null(object)) {
    stop("No object specified.")
  }

  if (is.null(exclude.assay)) {
    stop("No exon-excluded assay specified.")
  }

  exon.assay <- exon.assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- exon.assay
  
  df <- object[[exon.assay]][[]] 
  if (gene.name %ni% colnames(df)) {
    stop("No gene name found in meta data. Run ParseExonName for exon assay first.")
  }

  if (exon.name %ni% colnames(df)) {
    object[[exon.assay]][[exon.name]] <- rownames(object)
    df <- object[[exon.assay]][[]] 
  }
  
  if (!is.null(genes)) {
    genes0 <- CheckBindName(object, gene.name, assay = exon.assay)
    genes <- intersect(genes, genes0)
  }  else {
    genes <- CheckBindName(object, gene.name, assay = exon.assay)
  }

  df <- df[which(df[[gene.name]] %in% genes),]

  features <- df[['exon_name']]

  if (!is.null(ident.1) | !is.null(cells.1)) {
    tb <- AlternativeExpressionTest(object, ident.1 = ident.1, ident.2 = ident.2, cells.1 = cells.1, cells.2 = cells.2, assay = exon.assay, bind.assay = exclude.assay,
                                    bind.name = "exon_name", test.use = "PSI", min.pct = min.pct.exon, min.pct.bind.feature = min.pct, mode = 3, perm = perm, seed = seed, threads = threads)
    tb[[gene.name]] = df[tb$feature, gene.name]
  } else {
    tb <- FindAllAEMarkers(object, assay = exon.assay, bind.assay = exclude.assay, bind.name = "exon_name", test.use = "PSI", node = node, features = features,
                           return.thresh = return.thresh,
                           min.pct = min.pct.exon, min.pct.bind.feature = min.pct, mode = 3, perm = perm, seed = seed, threads = threads)
    tb[[gene.name]] = df[tb$feature, gene.name]
  }
  DefaultAssay(object) <- old.assay
  #tb <- tb[,-c("bind.feature","bind.feature.pct.1","bind.feature.pct.2")]
  tb
}
