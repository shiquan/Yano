#' @export
setMethod(f = "QuickRecipe0",
          signature = signature(counts = "SMatrix"),
          definition = function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                assay = NULL,
                                ...
                                ) {

            assay <- assay %||% "RNA"
            counts <- CreateSeuratObject(counts = counts, min.cells = min.cells,
                                         min.features = min.features, assay = assay)
            return(QuickRecipe0(counts, assay = assay, ...))
          })

#' @export
setMethod(f = "QuickRecipe0",
          signature = signature(counts = "Seurat"),
          definition = function(counts = NULL, scale.factor = 1e4,
                                assay = NULL,
                                ...
                                ) {
            assay <- assay %||% DefaultAssay(counts)
            message(paste0("Set default assay to ", assay))
            DefaultAssay(counts) <- assay
            counts <- NormalizeData(counts, normalization.method = "LogNormalize",
                                    scale.factor = scale.factor)
            counts
          })

ProcessDimReduc <- function(object = NULL, ndim=20, resolution = 0.5, nvar= 3000, features = NULL)
{
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nvar)
  
  features <- features %||% VariableFeatures(object)
  features <- intersect(features,rownames(object))
  
  object <- ScaleData(object, features = features)
  object <- RunPCA(object, features = features)
  object <- FindNeighbors(object, dims = 1:ndim)
  object <- FindClusters(object, resolution = resolution)
  object <- RunUMAP(object, dims = 1:ndim)
  object
}
#' Quick clust single cell gene expression matrix by Seurat pipeline
#'
#' @rdname QuickRecipe
#' @import Seurat
#' @import Matrix
#'
#' @export
setMethod(f = "QuickRecipe",
          signature = signature(counts = "Seurat"),
          definition = function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                nvar = 3000, resolution = 0.5, assay = NULL,
                                ndim = 20, ...
                                ) {
            
            object <- QuickRecipe0(counts=counts, meta.data = meta.data,
                                   min.cells =min.cells, min.features = min.features,
                                   assay = NULL, ...)

            ProcessDimReduc(object, ndim=ndim, resolution=resolution, nvar=nvar)
          })

setMethod(f = "QuickRecipe",
          signature = signature(counts = "SMatrix"),
          definition = function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                nvar = 3000, resolution = 0.5, assay = NULL,
                                ndim = 20, ...
                                ) {
            
            object <- QuickRecipe0(counts=counts, meta.data = meta.data,
                                   min.cells =min.cells, min.features = min.features,
                                   assay = assay, ...)

            ProcessDimReduc(object, ndim=ndim, resolution=resolution, nvar=nvar)
          })
#' Build cell-cell weight matrix on cell lineage.
#' @param order.cells Predefined cell ranks, used for cell lineage analysis.
#' @param k.nn k-nearest neighbors.
#' @param self.weight Diagnoal value in the weight matrix.
#' @param scale Scale the weight values by row.
#' @returns A sprase weight matrix.
#' @importFrom Matrix sparseMatrix
#' @export
WeightLineage <- function(order.cells = NULL, k.nn = 9, self.weight = 0, scale = FALSE) {
  if (is.null(order.cells)) stop("No order.cells specified.")
  
  n <- length(order.cells)
  i <- .Call("build_weight_i", n, k.nn)

  W <- sparseMatrix(i = i + 1, j = rep(1:n, each =k.nn), x=1, dims =c(n,n))

  W@x[W@x > 1] <- 1
  diag(W) <- self.weight

  if (scale) W <- W/rowSums(W)

  colnames(W) <- order.cells
  rownames(W) <- order.cells
  W
}
#' Calcualte cell-cell weight matrix by feature space (usually PCA) or spatial locations.
#' @param object Seurat object.
#' @param reduction Which dimesional reduction (PCA/ICA/harmony/spatial) to use for the weight matrix. Default is "pca".
#' @param dims Which dimensions to use as input features.
#' @param k.nn k-nearest neighbors.
#' @param spatial Use tissue coordinates instead of pca. Require input object is spatially resolved data.
#' @param self.weight Diagnoal value in the weight matrix.
#' @param scale Scale the weight values by row.
#' @param order.cells Predefined cell ranks, used for cell lineage analysis.
#' @param cells If set, run weight matrix on these cells only.
#' @returns A sparse weight matrix.
#' @export
GetWeights <- function(object= NULL,                       
                       reduction = "pca",
                       dims = NULL,
                       k.nn = 9,
                       spatial = FALSE,
                       self.weight = 0,
                       scale = FALSE,
                       order.cells = NULL,
                       cells = NULL)
{
  if (!is.null(order.cells)) {
    return(WeightLineage(order.cells = order.cells, k.nn = k.nn, self.weight=self.weight, scale=scale))
  }
  
  cells <- cells %||% colnames(object)
  
  if (isTRUE(spatial)) {
    message("Build weights on tissue coordiantes")
    emb <- GetTissueCoordinates(object)
  } else {
    message(paste0("Build weights on ", reduction))
    emb <- Embeddings(object,reduction = reduction)
    if (!is.null(dims)) emb <- emb[,dims]
  }

  emb <- emb[cells,]
  
  knn.rlt <- nabor::knn(data=emb, query = emb, k=k.nn)
  ncell <- length(cells)
  W <- sparseMatrix(i = rep(c(1:ncell), k.nn), 
                    j = c(knn.rlt$nn.idx),
                    x = 1,
                    dims = c(ncell, ncell))

  diag(W) <- self.weight
  
  if (scale) W <- W/rowSums(W)
  W[is.na(W)] <- 0
  colnames(W) <- cells
  rownames(W) <- cells
  W
}
#' Calculate Moran's index for features in parallel.
#' @param object
#' @param assay
#' @param slot
#' @param spatial
#' @param W
#' @param scale.weight
#' @param reduction
#' @param dims
#' @param k.nn
#' @param order.cells
#' @param cells
#' @param features
#' @param threads
#' @import Matrix
#' @export
RunAutoCorr <- function(object = NULL, 
                        assay = NULL,                        
                        slot = "data",
                        spatial = FALSE,
                        W = NULL,                              
                        scale.weight = FALSE,
                        reduction = "pca",
                        dims = NULL,
                        k.nn = 9,
                        order.cells = NULL,
                        cells = NULL,
                        features = NULL,
                        #perm = 10000,
                        threads=0)
{
  assay <- assay %||% DefaultAssay(object = object)
  message(paste0("Working on assay : ", assay))
  cells <- cells %||% colnames(object)
  cells <- intersect(colnames(object), cells)

  features <- features %||% rownames(object)
  features <- intersect(rownames(object),features)

  threads <- getCores(threads)
  
  if (is.null(W)) {
    W <- GetWeights(object=object, reduction=reduction,
                    dims=dims,
                    k.nn=k.nn,
                    cells=cells,
                    order.cells=order.cells,
                    self.weight = 0,
                    scale=TRUE,
                    spatial=spatial)
  } else {
    dims <- dim(W)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
    if (dims[1] != length(cells)) {
      cells <- intersect(colnames(W),cells)
      W <- W[cells, cells]
      diag(W) <- 0
      if (scale.weight) W <- W/rowSums(W)
    }
  }

  if (!is.null(order.cells)) {
    message("Working on lineage trajectory mode ..")
    cells <- order.cells    
  }
  x0 <- GetAssayData(object, assay = assay, slot = slot)[,cells]
  features <- intersect(rownames(x0), features)
  x0 <- x0[features,]
  x0 <- as(x0, "dgCMatrix")
  W <- as(W, "dgCMatrix")
  
  message(paste0("Run autocorrelation test for ", length(features), " features."))
  moransi.vals <- .Call("autocorrelation_test", x0, W, TRUE, threads);

  if (length(moransi.vals) == 1) stop(moransi.vals[[1]])
  
  #message("Run permutation test.")
  #moransi.vals2 <-  .Call("moransi_mc_test", x0, W, TRUE, perm, threads);

  ## message(paste0("Run Geary's C for ", length(features), " features."))
  ## gearysc.vals <- .Call("GearysC_test", x0, W);
  Ivals <- moransi.vals[[1]]
  IZvals <- moransi.vals[[3]]  
  names(Ivals) <- features
  pvals <- pnorm(IZvals, lower.tail = FALSE)
  names(pvals) <- features
  object[[assay]]@meta.features[['MoransI']] <- Ivals[rownames(object)]
  object[[assay]]@meta.features[['MoransI.pval']] <- pvals[rownames(object)]
  rm(W)
  rm(moransi.vals)
  gc()
  
  object
}

#' Set autocorrection by rank and/or score
#' @export
SetAutoCorrFeatures <- function(object = NULL,
                                moransi.min = 0,
                                assay = DefaultAssay(object),
                                p.cutoff = 0.05)
{
  tab <- object[[assay]]@meta.features

  cn <- colnames(tab)
  if ("MoransI" %ni% cn | "MoransI.pval" %ni% cn) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }

  idx <- which(tab[["MoransI"]] > moransi.min & tab[["MoransI.pval"]] <= p.cutoff)

  message(paste0(length(idx), " autocorrelated features."))
  tab[["AutoCorrFeature"]] <- FALSE
  tab[idx,][["AutoCorrFeature"]] <- TRUE
  object[[assay]]@meta.features <- tab
  
  object
}
#' @export
AutoCorrFeatures <- function(object = NULL, assay = NULL)
{
  assay <- assay %||% DefaultAssay(object)
  tab <- object[[assay]]@meta.features

  if ('AutoCorrFeature' %ni% colnames(tab)) {
    stop("No autocorrelation flag found, use SetAutoCorrFeature first.")
  }
  rownames(object)[which(tab[['AutoCorrFeature']])]
}


"%ni%" <- Negate("%in%")
#'
#' @import Matrix
#' @export
LocalCorr <- function(object = NULL,
                      moransi.cutoff = 0,
                      features = NULL,
                      cells = NULL,
                      assay = DefaultAssay(object),             
                      slot = "data",
                      scale = FALSE,
                      W = NULL,
                      reduction = "pca",
                      dims=NULL,
                      k.nn = 9,
                      clust.method = "ward.D2",
                      self.weight = 1,
                      verbose = TRUE
                      ) {

  tab <- object[[assay]]@meta.features
  
  if ("MoransI.value" %ni% colnames(tab) | "MoransI.rank" %ni% colnames(tab)) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }

  features <- features %||%  AutoCorrFeatures(object)
  features <- features %||%  rownames(object)  
  features <- intersect(features, rownames(tab))

  cells <- cells %||% colnames(object)
  cells <- intersect(cells, colnames(object))
  
  mtx <- GetAssayData(object, assay = assay, slot = slot)
  ## in case some features no coverage
  mtx <- mtx[which(rowSums(mtx) > 0),]

  features <- intersect(rownames(mtx), features)
  mtx <- mtx[features,]

  if (is.null(W)) {
    W <- GetWeights(object = object, reduction = reduction,
                    dims = dims,
                    k.nn = k.nn,
                    self.weight = self.weight,
                    scale = TRUE,
                    cells = cells)
  } else {
    dims <- dim(W)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
    if (dims[1] != length(cells)) {
      cells <- intersect(colnames(W),cells)
      W <- W[cells, cells]
      diag(W) <- 0
      if (scale) W <- W/rowSums(W)
    }
  }

  # mtx <- mtx %*% W
  mtx <- Matrix::tcrossprod(mtx,W)

  hm <- Matrix::tcrossprod(mtx)
  H <- hm/(length(cells)-1)
  H <- as.matrix(H)
  rm(mtx)
  rm(hm)
  gc()
  d <- dist(H)
  hc <- hclust(d, method = clust.method)
  H <- H[hc$order, hc$order]
  r <- SimpleList(LC = H, dist = d, hclust = hc)
  r
}
#'
#' @import RColorBrewer
#' @importFrom pheatmap pheatmap
#' @importFrom gtools mixedsort
#' @export
GroupLocalCorr <- function(lc = NULL, k = 10, plot = TRUE, name = "module")
{
  ## todo
  mod <- cutree(lc$hclust,k=k)
  mod.names <- names(mod)
  nm <- paste0(name,mod)
  nm <- factor(nm, levels = mixedsort(unique(nm)))
  ann <- data.frame(module=nm, row.names = mod.names)
  
  lc$module <- ann
  if (plot) {
    fig <- pheatmap(lc$LC,
                    cluster_rows = FALSE,
                    cluster_cols =FALSE,
                    show_rownames = FALSE,
                    show_colnames = FALSE,
                    annotation_row = ann)
    fig
  }
  lc
}
#'@importFrom Seurat AddModuleScore
#'@export
AddLCModule <- function(object = NULL, lc = NULL, min.features.per.module = 10, module.prefix.name = "module")
{
  #lc$module[["name"]] <- paste0(name, lc$module[,1])
  ml <- split(rownames(lc$module), lc$module[['module']])
  sel <- which(lapply(ml, length) >= min.features.per.module)
  ml <- ml[sel]

  meta <- setdiff(colnames(object@meta.data), unique(lc$module[['module']]))
  
  ## remove old modules
  object@meta.data <- object@meta.data[,meta]
  
  object <- AddModuleScore(object, features=ml, name=module.prefix.name)
  object
}

#'@importFrom parallel detectCores
getCores <- function(threads = 0)
{
  if (threads > 0) return(threads)
  return(detectCores())
}

#' @export
ParseExonName <- function(object = NULL, assay = NULL, stranded = TRUE)
{
  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay

  nm <- rownames(object)

  object[[assay]]@meta.features[['chr']] <- gsub("(.*):(.*)-(.*)/(.)/(.*)","\\1",nm)
  object[[assay]]@meta.features[['start']] <- as.integer(gsub("(.*):(.*)-(.*)/(.)/(.*)","\\2",nm))
  object[[assay]]@meta.features[['end']] <- as.integer(gsub("(.*):(.*)-(.*)/(.)/(.*)","\\3",nm))
  object[[assay]]@meta.features[['strand']] <- gsub("(.*):(.*)-(.*)/(.)/(.*)","\\4",nm)
  object[[assay]]@meta.features[['gene_name']] <- gsub("(.*):(.*)-(.*)/(.)/(.*)","\\5",nm)

  DefaultAssay(object) <- old.assay
  return(object)
}
#'
#' @importFrom data.table fread
#' @export
LoadEPTanno <- function(file = NULL, object = NULL, assay = NULL, stranded = TRUE)
{
  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  
  bed <- fread(file)[,c(1:9)]
  colnames(bed) <- c("chr","start","end","name","score","strand","n_gene","gene_name","type")
  
  bed$chr <- gsub("_","-",bed$chr)
  
  if (isTRUE(stranded)) {
    bed$name <- paste0(bed$chr,":",bed$start,"-",bed$end,"/",bed$strand)
  } else {
    bed$name <-  paste0(bed$chr,":",bed$start,"-",bed$end)
  }
  
  bed <- as.data.frame(bed)
  rownames(bed) <- bed$name

  features <- intersect(rownames(bed), rownames(object))

  if (length(features) == 0) {
    stop(paste0("No features found in assay ", assay))
  }

  message(paste0("Intersect ", length(features), " features."))
  
  bed <- bed[rownames(object),]
  
  object[[assay]]@meta.features[['chr']] <- bed[['chr']]
  object[[assay]]@meta.features[['start']] <- bed[['start']]
  object[[assay]]@meta.features[['end']] <- bed[['end']]
  object[[assay]]@meta.features[['name']] <- bed[['name']]
  object[[assay]]@meta.features[['strand']] <- bed[['strand']]
  object[[assay]]@meta.features[['n_gene']] <- bed[['n_gene']]
  object[[assay]]@meta.features[['gene_name']] <- bed[['gene_name']]
  object[[assay]]@meta.features[['type']] <- bed[['type']]

  DefaultAssay(object) <- old.assay

  object
}

#' @importFrom data.table fread 
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom stringr str_detect
#' @export
LoadVARanno <- function(file = NULL, object = NULL, assay = NULL, stranded = TRUE)
{
  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))

  bed <- fread(file)[,c(1:9)]
  colnames(bed) <- c("chr","start","end","name","score","strand","n_gene","gene_name","type")
  if (isTRUE(stranded)) {
    bed$name <- paste0(bed$chr,":",bed$start,"-",bed$end,"/",bed$strand)
  } else {
    bed$name <-  paste0(bed$chr,":",bed$start,"-",bed$end)
  }
  
  bed <- as.data.frame(bed)
  rownames(bed) <- bed$name

  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay


  if (isTRUE(stranded)) {
    locs <- gsub("(.*:[0-9]+)([ACGT=>]*).*/([-+])", "\\1/\\3",rownames(object))
    strands <- gsub(".*/([-+])","\\1",rownames(object))
  } else {
    locs <- gsub("(.*:[0-9]+)([ACGT=>]*).*", "\\1",rownames(object))
    strands <- "."
  }
  chrs <- gsub("(.*):.*","\\1",rownames(object))
  starts <- as.numeric(gsub("(.*):([0-9]+).*","\\2",rownames(object)))
    
  gv <- GRanges(chrs,IRanges(start=starts,width=1),strand=strands)
  gv$name <- rownames(object)

  gr <- GRanges(bed$chr,IRanges(start=bed$start,end=bed$end),strand=bed$strand)
  gr$name <- bed$name

  ov <- findOverlaps(gv, gr)
  var.sel <- gv[queryHits(ov)]$name
  ept.sel <- gr[subjectHits(ov)]$name
  names(ept.sel) <- var.sel
  gnames <- bed[ept.sel,]$gene_name
  types <- bed[ept.sel,]$type
  names(gnames) <- var.sel
  names(types) <- var.sel

  idx <- which(str_detect(rownames(object),"="))
  
  object[[assay]]@meta.features[['chr']] <- chrs
  object[[assay]]@meta.features[['start']] <- starts
  object[[assay]]@meta.features[['strand']] <- strands
  object[[assay]]@meta.features[['locus']] <- locs
  object[[assay]]@meta.features[['ept']] <- ept.sel[rownames(object)]
  object[[assay]]@meta.features[['gene_name']] <- gnames[rownames(object)]
  object[[assay]]@meta.features[['ept_type']] <- types[rownames(object)]
  object[[assay]]@meta.features[['type']] <- "alt"

  object[[assay]]@meta.features[['type']][idx] <- "ref"
  
  DefaultAssay(object) <- old.assay

  object
}


#' @importFrom Matrix sparseMatrix
#' @export
RunBlockCorr <- function(object = NULL,
                         bind.name = "gene_name",
                         features = NULL,
                         assay = NULL,
                         bind.assay = NULL,
                         bind.features = NULL,
                         prefix = NULL,
                         cells = NULL,
                         order.cells = NULL,
                         feature.types = NULL,                         
                         min.features.per.block = 2,
                         scale.factor = 1e4,
                         reduction = "pca",
                         sensitive.mode = FALSE,
                         spatial = FALSE,
                         W = NULL,
                         cell.size = NULL,
                         dims = NULL,
                         k.nn = 9,
                         perm=1000,
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
  
  tt <- Sys.time()
  cells <- cells %||% order.cells
  cells <- cells %||% colnames(object)
  cells <- intersect(cells, colnames(object))
  
  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))

  prefix <- prefix %||% bind.name
  
  features <- features %||% AutoCorrFeatures(object)
  features <- intersect(features, rownames(object))

  message(paste0("Working on ", length(features), " features."))

  threads <- getCores(threads)
  
  # Make weights
  if (is.null(W)) {
    W <- GetWeights(object = object,
                    reduction = reduction,
                    dims=dims,
                    k.nn = k.nn,
                    order.cells = order.cells,
                    self.weight = 1,
                    spatial=spatial,
                    scale=TRUE)
  } else {
    dims <- dim(W)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
  } 

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
  ## blocks <- unique(blocks)

  x <- GetAssayData(object, assay = assay, slot = "counts")
  
  # cell sizes
  cs <- cell.size %||% colSums(x)  
  cells <- intersect(cells,names(which(cs > 0)))
  W <- W[cells, cells]
  cs <- cs[cells]
  x <- x[,cells]
  ncell <- length(cells)
  
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
  
  message("Test dissimlarity of two processes ..")
  gc()
  ta <- .Call("D_test", x, y, W, perm, threads, idx, bidx, cs, scale.factor, sensitive.mode);
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
  #tab[[paste0(prefix, ".Lx")]] <- Lx[rownames(object)]
  #tab[[paste0(prefix, ".Ly")]] <- Ly[rownames(object)]
  tab[[paste0(prefix, ".pval")]] <- pval[rownames(object)]
  #tab[[paste0(prefix, ".padj")]] <- padj[rownames(object)]
  object[[assay]]@meta.features <- tab

  rm(ta)
  gc()

  tt <- Sys.time()-tt
  
  message(paste0("Runtime : ",format(tt)));
  object
}

#' @importFrom Matrix sparseMatrix
#' @export
RunCellCorr <- function(object = NULL,
                        meta.name = "nCount_RNA",
                        features = NULL,
                        assay = NULL,
                        prefix = NULL,
                        cells = NULL,
                        order.cells = NULL,
                        scale.factor = 1e4,
                        reduction = "pca",
                        sensitive.mode = FALSE,
                        spatial = FALSE,
                        W = NULL,
                        cell.size = NULL,
                        dims = NULL,
                        k.nn = 9,
                        perm=1000,
                        threads = 0
                        )
{
  cells <- cells %||% order.cells
  cells <- cells %||% colnames(object)
  cells <- intersect(cells, colnames(object))

  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  
  prefix <- prefix %||% meta.name
  
  features <- features %||% AutoCorrFeatures(object)
  features <- intersect(features, rownames(object))

  threads <- getCores(threads)
  
  message(paste0("Working on ", length(features), " features."))
  
  # Make weights
  if (is.null(W)) {
    W <- GetWeights(object = object,
                    reduction = reduction,
                    dims=dims,
                    k.nn = k.nn,
                    order.cells = order.cells,
                    self.weight = 1,
                    spatial=spatial,
                    scale=TRUE)
  } else {
    dims <- dim(W)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
  } 

  tab <- object@meta.data

  if (meta.name %ni% colnames(tab)) {
    stop(paste0("No meta.name found in the meta.data."))
  }
  
  tab <- tab[which(tab[[meta.name]] != "."),] # skip empty
  y <- tab[[meta.name]]
  names(y) <- rownames(tab)
  
  if (!is.numeric(y)) stop("meta.name contain non-numeric value.")
  
  features <- intersect(features, rownames(object))
    
  if (length(features) == 0) {
    stop("No features found.")
  }

  message(paste0("Processing ", length(cells), " cells.."))
  tab <- tab[cells,]
  
  x <- GetAssayData(object, assay = assay, slot = "counts")
  
  #cells <- intersect(cells,names(which(cs > 0)))
  W <- W[cells, cells]
  cs <- cs[cells]
  x <- x[,cells]
  y <- y[cells]
  ncell <- length(cells)
  
  idx <- match(features, rownames(x))

  cs <- cell.size %||% colSums(x)
  message("Test dissimlarity of two processes ..")
  gc()
  ta <- .Call("D_test_cell", x, y, W, perm, threads, idx, cs, scale.factor, sensitive.mode);
  if (length(ta) == 1) stop(ta[[1]])

  Lx <- ta[[1]]
  r <- ta[[2]]
  e <- ta[[3]]
  tval <- ta[[4]]

  pval <- pt(tval, df = perm - 1, lower.tail = FALSE)
  
  names(pval) <- features  
  names(Lx) <- features
  names(r) <- features
  names(e) <- features
  
  tab <- object[[assay]]@meta.features
  tab[[paste0(prefix, ".E")]] <- e[rownames(object)]
  tab[[paste0(prefix, ".r")]] <- r[rownames(object)]
  #tab[[paste0(prefix, ".Lx")]] <- Lx[rownames(object)]
  tab[[paste0(prefix, ".pval")]] <- pval[rownames(object)]
  object[[assay]]@meta.features <- tab

  rm(ta)
  gc()

  DefaultAssay(object) <- old.assay
  object
}

#'
#' @export
RunTwoAssayCorr <- function(object = NULL,
                            assay1 = "RNA",
                            assay2 = "EPT",
                            features1 = VariableFeatures(object, assay = RNA),
                            features2 = AutoCorrFeatures(object, assay = assay2),
                            cells = NULL,
                            slot = "data",
                            W = NULL,
                            scale = FALSE,
                            reduction = "pca",
                            dims=NULL,
                            k.nn = 5,
                            self.weight = 1,
                            verbose = TRUE
                            )
{
  cells <- cells %||% colnames(object)
  cells <- intersect(colnames(object), cells)

  if (is.null(features1)) {
    stop("Features1 is empty.")
  }

  if (is.null(features2)) {
    stop("Features2 is empty.")
  }

  # Make weights
  if (is.null(W)) {
    W <- GetWeights(object = object, reduction = reduction, dims=dims, k.nn = k.nn,
                    self.weight = self.weight,
                    cells = cells)
  } else {
    dims <- dim(W)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
    if (dims[1] != length(cells)) {
      cells <- intersect(colnames(W),cells)
      W <- W[cells, cells]
      diag(W) <- 0
      if (scale) W <- W/rowSums(W)
    }
  }  
  
  old.assay <- DefaultAssay(obj)
  
  DefaultAssay(obj) <- assay1  
  features1 <- intersect(rownames(object),features1)
  if (is.null(features1)) {
    stop("No valid features1 found, check the assay1 name.")
  }
  x <- FetchData(obj, vars = features1, cells = cells, slot = slot)
  
  DefaultAssay(obj) <- assay2
  features2 <- intersect(rownames(object),features2)
  if (is.null(features2)) {
    stop("No valid features2 found, check the assay2 name.")
  }
  y <- FetchData(obj, vars = features2, cells = cells, slot = slot)

  x <- tcrossprod(x, W)
  y <- tcrossprod(y, W)
  x <- t(scale(t(x)))
  y <- t(scale(t(y)))
  
  # H(x) = sum(X,Y)/sqrt(sum(X^2)*sum(Y^2))  
  hm <- x %*% t(y)

  H <- hm/(length(cells)-1)
  H <- as.matrix(H)
  rm(x)
  rm(y)
  rm(hm)
  gc()
  d <- dist(H)
  hc <- hclust(d, method = method)
  H <- H[hc$order, hc$order]
  r <- SimpleList(LC = H, dist = d, hclust = hc)
  r
}

#' @export
aggregateGeneSets <- function(object = NULL, name = NULL, features = NULL, assay = NULL, scale.factor=1e6, scaled=TRUE)
{
  if (is.null(name)) {
    stop("Set up name for the set.")
  }

  if (is.null(features)) {
    stop("Set up features for the set.")
  }

  assay <- assay %||%  DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  features <- intersect(features, rownames(object)) 

  slot <- "counts"

  dat <- GetAssayData(object, slot=slot)
  cs  <- colSums(dat)
  dat <- dat[features,]
  dat <- colSums(dat)
  
  dat <- log1p(dat/cs*scale.factor)

  if (scaled) dat <- scale(dat)
  object[[name]] <- dat

  DefaultAssay(object) <- old.assay
  
  object
}


#' @export
aggregateCellByGroup <- function(object = NULL, cell.group = NULL, features = NULL, assay = NULL, avg.by.cells = FALSE)
{
  cell.group <- cell.group %||% Idents(object)
  assay <- assay %||%  DefaultAssay(object)
  DefaultAssay(object) <- assay
  features <- features %||% rownames(object)
  features <- intersect(features, rownames(object)) 

  slot <- "counts"

  dat <- GetAssayData(object, slot=slot)
  dat <- dat[features,]
  groups <- unique(cell.group)
  dat <- as(dat, "dgTMatrix")

  m <- Matrix::sparseMatrix(
    i = dat@i+1,
    j = match(cell.group[dat@j+1], groups),
    x = dat@x,
    dims = c(length(features), length(groups)))
  
  colnames(m) <- groups
  rownames(m) <- features

  if (avg.by.cells) {
    cs <- colSums(m)
    m <- m/cs
  }
  
  m
}
#'@importFrom dplyr %>%
#'@importFrom gtools mixedsort
#'@importFrom ggrepel geom_label_repel
#'@import ggplot2
#'@import ggrepel
#'@importFrom scattermore geom_scattermore
#' 
#'@export
FbtPlot <- function(object = NULL, assay = NULL, chr = "chr", start = "start", val = NULL, col.by = NULL, cols = NULL, sel.chrs = NULL, xlab = "Chromosome", ylab = expression(-log[10](P)), raster = NULL, types = NULL, point.label = NULL, arrange.type = FALSE, label.size=3,...)
{
  assay <- assay %||% DefaultAssay(object)
  tab0 <- object[[assay]]@meta.features
  cols <- cols %||% c("#131313","blue","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
  
  if (is.null(val)) stop("No value name specified.")  
  if (chr %ni% colnames(tab0)) stop("No chr name found.")
  if (start %ni% colnames(tab0)) stop("No start name found.")
  if (val %ni% colnames(tab0)) stop("No val name found.")
  
  if (!is.null(sel.chrs)) {
    tab0 <- tab0 %>% filter (chr %in% sel.chrs)
  }

  if (!is.null(types)) {
    if ("type" %ni% colnames(tab0)) stop("No type found.")
    tab0 <- subset(tab0, type %in% types)
    if (nrow(tab0) == 0) stop("Empty records.")
  }
  lv <- mixedsort(unique(tab0[[chr]]))
  
  tab <- data.frame(chr = factor(tab0[[chr]], levels = lv),
                    start = as.numeric(tab0[[start]]),
                    pval = -log10(as.numeric(tab0[[val]])),
                    row.names = rownames(tab0))
  
  if (!is.null(col.by)) {
    tab[[col.by]] <- tab0[[col.by]]
  }

  tab <- subset(tab, !is.na(pval))
  
  data_cum <- tab %>% group_by(chr) %>% summarise(max_bp = max(start)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% select(chr, bp_add)
  
  data <- tab %>% inner_join(data_cum, by = "chr") %>% mutate(bp_cum = start + bp_add)
  axis_set <- data %>% group_by(chr) %>% summarize(center = mean(bp_cum))

  data$name <- rownames(tab)
  rownames(data) <- rownames(tab)
  
  if (isTRUE(arrange.type)) data <- data %>% arrange(type)

  fbt_theme <- function() {
    theme(
      legend.text = element_text(face = "italic",color = "black",family = "Helvetica",size = rel(1.5)),
      axis.title.y = element_text(color = "black", family = "Helvetica",size = rel(1)),
      axis.title.x = element_text(color = "black", family = "Helvetica",size = rel(1.5)),
      axis.text = element_text(family = "Helvetica",color = "black",size = rel(1.5)),
      axis.line = element_blank(),
      axis.ticks = element_line(size = rel(1), color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size= rel(2), linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "whitesmoke"),
      legend.key = element_rect(fill = "whitesmoke"),
      legend.title = element_text(size = rel(1.5),family = "Helvetica"),
      plot.title = element_text(color = "black",face = "bold",size = rel(1.7),family = "Helvetica")
    )
  }

  xi <- data_cum$bp_add
  xi <- xi[-1]
  p <- ggplot(data) + geom_vline(xintercept = xi, color="red", linetype="dotted")

  raster <- raster %||% (nrow(x = data) > 1e5)
  
  if (!is.null(col.by)) {
    if (raster) {
      p <- p + geom_scattermore(mapping = aes(x=bp_cum, y = pval, col = .data[[col.by]]),shape=21,...)
      p <- p + scale_color_manual(values = cols)
    } else {
      p <- p + geom_point(aes(x=bp_cum, y=pval, fill=.data[[col.by]]), shape=21,...)
      p <- p + scale_fill_manual(values = cols)
    }


  } else {
    if (raster) {
      p <- p + geom_scattermore(aes(x=bp_cum, y=pval),  ...)    
    } else {
      p <- p + geom_point(aes(x=bp_cum, y=pval),  ...)    
    }
  }
  p <- p + scale_x_continuous(label = axis_set$chr, breaks = axis_set$center,
                              limits = c(min(data$bp_cum), max(data$bp_cum)),
                              expand=c(0,0))
  p <- p + fbt_theme() + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))
  p <- p + xlab(xlab) + ylab(ylab)

  if (!is.null(point.label)) {
    sel <- intersect(point.label, rownames(data))
    if (length(sel) > 0) {
      p <- p + geom_label_repel(data=data[sel,],aes(x=bp_cum, y=pval,label=name),box.padding = 0.5, max.overlaps = Inf, size=label.size)
    }
  }
  p
}
