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


"%ni%" <- Negate("%in%")
## #'
## #' @import Matrix
## #' @export
## LocalCorr <- function(object = NULL,
##                       moransi.cutoff = 0,
##                       features = NULL,
##                       cells = NULL,
##                       assay = DefaultAssay(object),             
##                       slot = "data",
##                       scale = FALSE,
##                       W = NULL,
##                       reduction = "pca",
##                       dims=NULL,
##                       k.nn = 9,
##                       clust.method = "ward.D2",
##                       self.weight = 1,
##                       verbose = TRUE
##                       ) {

##   tab <- object[[assay]]@meta.features
  
##   if ("MoransI.value" %ni% colnames(tab) | "MoransI.rank" %ni% colnames(tab)) {
##     stop("No Morans'I value found, use RunAutoCorr first.")
##   }

##   features <- features %||%  AutoCorrFeatures(object)
##   features <- features %||%  rownames(object)  
##   features <- intersect(features, rownames(tab))

##   cells <- cells %||% colnames(object)
##   cells <- intersect(cells, colnames(object))
  
##   mtx <- GetAssayData(object, assay = assay, slot = slot)
##   ## in case some features no coverage
##   mtx <- mtx[which(rowSums(mtx) > 0),]

##   features <- intersect(rownames(mtx), features)
##   mtx <- mtx[features,]

##   if (is.null(W)) {
##     W <- GetWeights(object = object, reduction = reduction,
##                     dims = dims,
##                     k.nn = k.nn,
##                     self.weight = self.weight,
##                     scale = TRUE)
##     #cells = cells)
##   } else {
##     dims <- dim(W)
##     if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
##     if (dims[1] != length(cells)) {
##       cells <- intersect(colnames(W),cells)
##       W <- W[cells, cells]
##       diag(W) <- 0
##       if (scale) W <- W/rowSums(W)
##     }
##   }

##   # mtx <- mtx %*% W
##   mtx <- Matrix::tcrossprod(mtx,W)

##   hm <- Matrix::tcrossprod(mtx)
##   H <- hm/(length(cells)-1)
##   H <- as.matrix(H)
##   rm(mtx)
##   rm(hm)
##   gc()
##   d <- dist(H)
##   hc <- hclust(d, method = clust.method)
##   H <- H[hc$order, hc$order]
##   r <- SimpleList(LC = H, dist = d, hclust = hc)
##   r
## }
## #'
## #' @import RColorBrewer
## #' @importFrom pheatmap pheatmap
## #' @importFrom gtools mixedsort
## #' @export
## GroupLocalCorr <- function(lc = NULL, k = 10, plot = TRUE, name = "module")
## {
##   ## todo
##   mod <- cutree(lc$hclust,k=k)
##   mod.names <- names(mod)
##   nm <- paste0(name,mod)
##   nm <- factor(nm, levels = mixedsort(unique(nm)))
##   ann <- data.frame(module=nm, row.names = mod.names)
  
##   lc$module <- ann
##   if (plot) {
##     fig <- pheatmap(lc$LC,
##                     cluster_rows = FALSE,
##                     cluster_cols =FALSE,
##                     show_rownames = FALSE,
##                     show_colnames = FALSE,
##                     annotation_row = ann)
##     fig
##   }
##   lc
## }
## #'@importFrom Seurat AddModuleScore
## #'@export
## AddLCModule <- function(object = NULL, lc = NULL, min.features.per.module = 10, module.prefix.name = "module")
## {
##   #lc$module[["name"]] <- paste0(name, lc$module[,1])
##   ml <- split(rownames(lc$module), lc$module[['module']])
##   sel <- which(lapply(ml, length) >= min.features.per.module)
##   ml <- ml[sel]

##   meta <- setdiff(colnames(object@meta.data), unique(lc$module[['module']]))
  
##   ## remove old modules
##   object@meta.data <- object@meta.data[,meta]
  
##   object <- AddModuleScore(object, features=ml, name=module.prefix.name)
##   object
## }

#'@importFrom parallel detectCores
getCores <- function(threads = 0)
{
  if (threads > 0) return(threads)
  threads <- detectCores() - 1
  if (threads > 1) return(threads)
  return(1)
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

## #' @importFrom Matrix sparseMatrix
## #' @export
## RunCellCorr <- function(object = NULL,
##                         meta.name = "nCount_RNA",
##                         features = NULL,
##                         assay = NULL,
##                         prefix = NULL,
##                         cells = NULL,
##                         order.cells = NULL,
##                         scale.factor = 1e4,
##                         reduction = "pca",
##                         sensitive.mode = FALSE,
##                         spatial = FALSE,
##                         W = NULL,
##                         cell.size = NULL,
##                         dims = NULL,
##                         k.nn = 9,
##                         perm=1000,
##                         threads = 0
##                         )
## {
##   cells <- cells %||% order.cells
##   cells <- cells %||% colnames(object)
##   cells <- intersect(cells, colnames(object))

##   assay <- assay %||% DefaultAssay(object)
##   message(paste0("Working on assay ", assay))
##   old.assay <- DefaultAssay(object)
##   DefaultAssay(object) <- assay
  
##   prefix <- prefix %||% meta.name
  
##   features <- features %||% AutoCorrFeatures(object)
##   features <- intersect(features, rownames(object))

##   threads <- getCores(threads)
  
##   message(paste0("Working on ", length(features), " features."))
  
##   # Make weights
##   if (is.null(W)) {
##     W <- GetWeights(object = object,
##                     reduction = reduction,
##                     dims=dims,
##                     k.nn = k.nn,
##                     order.cells = order.cells,
##                     self.weight = 1,
##                     spatial=spatial,
##                     scale=TRUE)
##   } else {
##     dims <- dim(W)
##     if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
##   } 

##   tab <- object@meta.data

##   if (meta.name %ni% colnames(tab)) {
##     stop(paste0("No meta.name found in the meta.data."))
##   }
  
##   tab <- tab[which(tab[[meta.name]] != "."),] # skip empty
##   y <- tab[[meta.name]]
##   names(y) <- rownames(tab)
  
##   if (!is.numeric(y)) stop("meta.name contain non-numeric value.")
  
##   features <- intersect(features, rownames(object))
    
##   if (length(features) == 0) {
##     stop("No features found.")
##   }

##   message(paste0("Processing ", length(cells), " cells.."))
##   tab <- tab[cells,]
  
##   x <- GetAssayData(object, assay = assay, slot = "counts")
  
##   #cells <- intersect(cells,names(which(cs > 0)))
##   W <- W[cells, cells]
##   cs <- cs[cells]
##   x <- x[,cells]
##   y <- y[cells]
##   ncell <- length(cells)
  
##   idx <- match(features, rownames(x))

##   cs <- cell.size %||% colSums(x)
##   message("Test dissimlarity of two processes ..")
##   gc()
##   ta <- .Call("D_test_cell", x, y, W, perm, threads, idx, cs, scale.factor, sensitive.mode);
##   if (length(ta) == 1) stop(ta[[1]])

##   Lx <- ta[[1]]
##   r <- ta[[2]]
##   e <- ta[[3]]
##   tval <- ta[[4]]

##   pval <- pt(tval, df = perm - 1, lower.tail = FALSE)
  
##   names(pval) <- features  
##   names(Lx) <- features
##   names(r) <- features
##   names(e) <- features
  
##   tab <- object[[assay]]@meta.features
##   tab[[paste0(prefix, ".E")]] <- e[rownames(object)]
##   tab[[paste0(prefix, ".r")]] <- r[rownames(object)]
##   #tab[[paste0(prefix, ".Lx")]] <- Lx[rownames(object)]
##   tab[[paste0(prefix, ".pval")]] <- pval[rownames(object)]
##   object[[assay]]@meta.features <- tab

##   rm(ta)
##   gc()

##   DefaultAssay(object) <- old.assay
##   object
## }

## #'
## #' @export
## RunTwoAssayCorr <- function(object = NULL,
##                             assay1 = "RNA",
##                             assay2 = "EPT",
##                             features1 = VariableFeatures(object, assay = RNA),
##                             features2 = AutoCorrFeatures(object, assay = assay2),
##                             cells = NULL,
##                             slot = "data",
##                             W = NULL,
##                             scale = FALSE,
##                             reduction = "pca",
##                             dims=NULL,
##                             k.nn = 5,
##                             self.weight = 1,
##                             verbose = TRUE
##                             )
## {
##   cells <- cells %||% colnames(object)
##   cells <- intersect(colnames(object), cells)

##   if (is.null(features1)) {
##     stop("Features1 is empty.")
##   }

##   if (is.null(features2)) {
##     stop("Features2 is empty.")
##   }

##   # Make weights
##   if (is.null(W)) {
##     W <- GetWeights(object = object, reduction = reduction, dims=dims, k.nn = k.nn,
##                     self.weight = self.weight)
##     #cells = cells)
##   } else {
##     dims <- dim(W)
##     if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
##     if (dims[1] != length(cells)) {
##       cells <- intersect(colnames(W),cells)
##       W <- W[cells, cells]
##       diag(W) <- 0
##       if (scale) W <- W/rowSums(W)
##     }
##   }  
  
##   old.assay <- DefaultAssay(obj)
  
##   DefaultAssay(obj) <- assay1  
##   features1 <- intersect(rownames(object),features1)
##   if (is.null(features1)) {
##     stop("No valid features1 found, check the assay1 name.")
##   }
##   x <- FetchData(obj, vars = features1, cells = cells, slot = slot)
  
##   DefaultAssay(obj) <- assay2
##   features2 <- intersect(rownames(object),features2)
##   if (is.null(features2)) {
##     stop("No valid features2 found, check the assay2 name.")
##   }
##   y <- FetchData(obj, vars = features2, cells = cells, slot = slot)

##   x <- tcrossprod(x, W)
##   y <- tcrossprod(y, W)
##   x <- t(scale(t(x)))
##   y <- t(scale(t(y)))
  
##   # H(x) = sum(X,Y)/sqrt(sum(X^2)*sum(Y^2))  
##   hm <- x %*% t(y)

##   H <- hm/(length(cells)-1)
##   H <- as.matrix(H)
##   rm(x)
##   rm(y)
##   rm(hm)
##   gc()
##   d <- dist(H)
##   hc <- hclust(d, method = method)
##   H <- H[hc$order, hc$order]
##   r <- SimpleList(LC = H, dist = d, hclust = hc)
##   r
## }

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
