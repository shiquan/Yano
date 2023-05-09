#' @export
setMethod(f = "QuickRecipe0",
          signature = signature(counts = "SMatrix"),
          definition = function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                nvar = 2000, assay = "RNA",
                                ...
                                ) {            
            counts <- CreateSeuratObject(counts = counts, min.cells = min.cells,
                                         min.features = min.features, assay = assay)
            return(QuickRecipe0(counts, assay = assay, nvar=nvar, ...))
          })

#' @export
setMethod(f = "QuickRecipe0",
          signature = signature(counts = "Seurat"),
          definition = function(counts = NULL, nvar = 2000, assay = NULL, scale.factor = 1e4,
                                ...
                                ) {
            assay <- assay %||% DefaultAssay(counts)
            message(paste0("Set default assay to ", assay))
            DefaultAssay(counts) <- assay

            counts <- NormalizeData(counts, normalization.method = "LogNormalize",
                                    scale.factor = scale.factor)
            
            counts <- FindVariableFeatures(counts, selection.method = "vst", nfeatures = nvar)
            counts
          })
#'
#' @export
CreateEPAssay <- function(counts = NULL,
                          assay = 'EP',
                          meta.data = NULL,
                          min.cells = 100) 
{
  counts <- CreateAssayObject(counts = counts, assay = assay, min.cells = min.cells)

  message(paste0("Set default assay to ", assay))
  DefaultAssay(counts) <- assay

  counts
}
#'
#' @export
ProcessDimReduc <- function(object = NULL, ndim=20, resolution = 0.5, features = NULL)
{
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
                                nvar = 3000, resolution = 0.5, assay = "RNA",
                                ndim = 20, ...
                                ) {
            
            object <- QuickRecipe0(counts=counts, meta.data = meta.data,
                                   min.cells =min.cells, min.features = min.features,
                                   nvar = nvar, assay = assay, ...)

            ProcessDimReduc(object, ndim=ndim, resolution=resolution)
          })

setMethod(f = "QuickRecipe",
          signature = signature(counts = "SMatrix"),
          definition = function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                nvar = 3000, resolution = 0.5, assay = "RNA",
                                ndim = 20, ...
                                ) {
            
            object <- QuickRecipe0(counts=counts, meta.data = meta.data,
                                   min.cells =min.cells, min.features = min.features,
                                   nvar = nvar, assay = assay, ...)

            ProcessDimReduc(object, ndim=ndim, resolution=resolution)
          })

#'
#' @export
GetWeights <- function(object= NULL,
                       reduction = "pca",
                       dims = NULL,
                       k.nn = 20,
                       kernel.method = "dist",
                       cells = NULL)
{
  emb <- Embeddings(object,reduction = reduction)
  if (!is.null(dims)) emb <- emb[,dims]

  knn.rlt <- nabor::knn(data=emb, query = emb, k=k.nn)
  if (kernel.method == "average") x = 1
  else x = 1/c(knn.rlt$nn.dists)

  if (is.null(cells)) cells <- colnames(object)
  ncell <- length(cells)
  W <- sparseMatrix(i = rep(c(1:ncell), k.nn), 
                    j = c(knn.rlt$nn.idx),
                    x = x,
                    dims = c(ncell, ncell))
  if (kernel.method == "dist") diag(W) <- 0
  else diag(W) <- 1
  
  W <- W/rowSums(W)
  W
}

#'
#' @import sparseMatrixStats
#' @import Matrix
#' @export
RunAutoCorr <- function(object = NULL, 
                        assay = DefaultAssay(object = object),
                        slot = "data",
                        scaled = FALSE,
                        weights = NULL,                              
                        weights.scaled = FALSE,
                        reduction = "pca",
                        dims = NULL,
                        k.nn = 20,
                        kernel.method = "dist",
                        cells = NULL,
                        features = NULL,
                        verbose = TRUE)
{
  cells <- cells %||% colnames(object)
  cells <- intersect(colnames(object), cells)

  features <- features %||% rownames(object)
  features <- intersect(rownames(object),features)
  
  if (is.null(weights)) {
    W <- GetWeights(object=object, reduction=reduction,dims=dims,k.nn=k.nn,kernel.method=kernel.method,
                    cells=cells)
  } else {
    dims <- dim(weights)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
    if (dims[1] != length(cells)) {
      cells <- intersect(colnames(weights),cells)
      W <- weights[cells, cells]
      diag(W) <- 0
      if (!weights.scaled) W <- W/rowSums(W)
    }
  }
  
  message(paste0("Working on assay : ", assay))
  x <- GetAssayData(object, assay = assay, slot = slot)[,cells]

  # in case some features missing in scaled matrix
  features <- intersect(rownames(x), features)
  message(paste0("Run autocorrelation for ", length(features), " features."))
  x <- x[features,]
  
  coverage <- rowSums(x > 0)/length(cells)
  names(coverage) <- features
  
  if (!scaled) {
    x <- t(scale(t(x)))
  }
  
  y <- Matrix::t(Matrix::tcrossprod(W, x))
  y <- as.matrix(y)

  #todo
  corr <- lineup::corbetw2mat(t(y),t(x))
  vals <- corr * rowSds(y)/rowSds(x)
  rm(x)
  rm(y)
  
  vals <- sort(vals,decreasing = TRUE)
  tab <- data.frame(MoransI.value=vals,
                    MoransI.rank=1:length(vals),
                    coverage = coverage[names(vals)],
                    row.names=names(vals))
  rm(vals)
  gc()
  
  tab0 <- object[[assay]]@meta.features
  tab <- tab[rownames(object),]
  # remove NAs in rownames
  rownames(tab) <- rownames(object) 

  # Set autocorrection features for downstream analysis, here roughtly set Moran's I greater than 0
  tab[['AutoCorrFeature']] = FALSE
  idx <- which(tab[["MoransI.value"]] > 0)
  tab[idx,][["AutoCorrFeature"]] <- TRUE

  object[[assay]]@meta.features <- cbind(tab0,tab)
  object
}

#' Set autocorrection by rank and/or score
#' @export
SetAutoCorrFeatures <- function(object = NULL,
                                moransi.min = 0,
                                top.n = 500,
                                assay = DefaultAssay(object),
                                plot = TRUE
                                )
{
  tab <- object[[assay]]@meta.features

  if ("MoransI.value" %ni% colnames(tab) | "MoransI.rank" %ni% colnames(tab)) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }

  idx <- which(tab[["MoransI.value"]] > moransi.min & tab[["MoransI.rank"]] <= top.n)
  
  tab[["AutoCorrFeature"]] <- FALSE
  tab[idx,][["AutoCorrFeature"]] <- TRUE
  object[[assay]]@meta.features <- tab

  if (plot) {
    v <- min(tab[idx,][["MoransI.value"]])
    plot(density(tab[['MoransI.value']], na.rm = TRUE))
    abline(v=v, col="red")
  }
  
  object
}
#' @export
AutoCorrFeatures <- function(object = NULL, assay = DefaultAssay(object))
{
  tab <- object[[assay]]@meta.features

  if ("MoransI.value" %ni% colnames(tab) | "MoransI.rank" %ni% colnames(tab)) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }
  rownames(object)[which(tab[['AutoCorrFeature']])]
}


"%ni%" <- Negate("%in%")
#'
#' @import Matrix
#' @export
LocalCorr <- function(object = NULL,
                      moransi.cutoff = 0,
                      features = AutoCorrFeatures(object),
                      cells = NULL,
                      assay = DefaultAssay(object),             
                      slot = "data",
                      scaled = FALSE,
                      weights = NULL,
                      weights.scaled = FALSE,
                      reduction = "pca",
                      dims=NULL,
                      k.nn = 20,
                      kernel.method = "dist",                      
                      method="average",
                      verbose = TRUE
                      ) {

  tab <- object[[assay]]@meta.features
  
  if ("MoransI.value" %ni% colnames(tab) | "MoransI.rank" %ni% colnames(tab)) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }

  features <- features %||% rownames(tab)
  features <- intersect(features, rownames(tab))

  cells <- cells %||% colnames(object)
  cells <- intersect(cells, colnames(object))
  
  mtx <- GetAssayData(object, assay = assay, slot = slot)
  features <- intersect(rownames(mtx), features)
  mtx <- mtx[features,]

  if (is.null(weights)) {
    W <- GetWeights(object = object, reduction = reduction,
                    dims = dims,
                    k.nn = k.nn,
                    kernel.method = kernel.method,
                    cells = cells)
  } else {
    dims <- dim(weights)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
    if (dims[1] != length(cells)) {
      cells <- intersect(colnames(weights),cells)
      W <- weights[cells, cells]
      diag(W) <- 0
      if (!weights.scaled) W <- W/rowSums(W)
    }
  }

  mtx <- mtx %*% W

  # smoothed by weights for logcounts, then scaled 
  if (!scaled) {
    mtx <- t(scale(t(mtx)))
  }

  # H(x) = sum(X,Y)/sqrt(sum(X^2)*sum(Y^2))  
  hm <- Matrix::tcrossprod(mtx) #mtx.2 %*% t(mtx.2)
  #rs <- rowSums(mtx * mtx)
  #rs2 <- sqrt(as.matrix(rs) %*% rs)  
  H <- hm/(length(cells)-1)
  H <- as.matrix(H)
  rm(mtx)
  rm(hm)
  gc()
  d <- dist(H)
  hc <- hclust(d, method = method)
  H <- H[hc$order, hc$order]
  r <- SimpleList(LC = H, dist = d, hclust = hc)
  r
}
#'
#' @import RColorBrewer
#' @import pheatmap
#' @export
GroupLocalCorr <- function(lc = NULL, k = 10, plot = TRUE, name="module")
{
  ## todo
  mod <- cutree(lc$hclust,k=k)
  mod.names <- names(mod)
  groups <- paste0(name,mod)

  ann <- data.frame(module=groups, row.names = mod.names)
  lc$module <- ann
  if (plot) {
    # require(RColorBrewer)
    # require(pheatmap)
    ## if (is.null(annotation_colors)) {
    ##   annotation_colors <- brewer.pal(12, "Paired")
    ## }
    ## if (k > length(annotation_colors)) {
    ##   pals <- colorRampPalette(annotation_colors)
    ##   annotation_colors <- pals(k)
    ## }

    fig <- pheatmap(lc$LC,
                    cluster_rows = FALSE,
                    cluster_cols =FALSE,
                    show_rownames = FALSE,
                    show_colnames = FALSE,
                    annotation_row = ann)
                    #annotation_colors = list(module=annotation_colors))
    fig
  }
  lc
}
#'@export
AddLCModule <- function(object = NULL, lc = NULL, min.features.per.module = 10)
{
  #lc$module[["name"]] <- paste0(name, lc$module[,1])
  ml <- split(rownames(lc$module), lc$module[['module']])
  sel <- which(lapply(ml, length) >= min.features.per.module)
  ml <- ml[sel]

  meta <- setdiff(colnames(object@meta.data), unique(lc$module[['module']]))
  
  ## remove old modules
  object@meta.data <- object@meta.data[,meta]
  
  object <- AddModuleScore(object, features=ml, name=name)
  object
}
#'
#' @import data.table
#' @export
LoadEPAnno <- function(file = NULL, object = NULL, assay = DefaultAssay(object))
{
  bed <- fread(file)[,c(1:9)]
  colnames(bed) <- c("chr","start","end","name","score","strand","n_gene","gene_name","type")
  bed$name <- paste(bed$chr,bed$start,bed$end,bed$strand,sep="-")
  bed <- as.data.frame(bed)
  rownames(bed) <- bed$name

  features <- intersect(rownames(bed), rownames(object))

  if (length(features) == 0) {
    stop(paste0("No features found in assay ", assay))
  }
  
  bed <- bed[rownames(object),]
  tab <- object[[assay]]@meta.features  
  nm <- setdiff(colnames(tab), colnames(bed))
  tab <- tab[,nm]
  
  object[[assay]]@meta.features <- cbind(tab,bed)
  object
}

#'
#' @export
RunEPBlockCorr <- function(object = NULL,
                           block.name = "gene_name",
                           name = "altDisp",
                           assay = "EP",
                           slot = "data",
                           features = AutoCorrFeatures(object),
                           alter.splice.mode = FALSE,
                           block.assay = "RNA",
                           block.force.replace = FALSE,
                           cells = NULL,                           
                           min.features.per.block = 2,
                           scale.factor = 1e4,
                           weights = NULL,
                           weights.scaled = FALSE,
                           reduction = "pca",
                           dims=NULL,
                           k.nn = 20,
                           kernel.method = "dist",
                           verbose = TRUE
                           )
{
  cells <- cells %||% colnames(object)
  cells <- intersect(colnames(object), cells)

  # Make weights
  if (is.null(weights)) {
    W <- GetWeights(object = object, reduction = reduction, dims=dims, k.nn = k.nn,
                    kernel.method = kernel.method,
                    cells = cells)
  } else {
    dims <- dim(weights)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
    if (dims[1] != length(cells)) {
      cells <- intersect(colnames(weights),cells)
      W <- weights[cells, cells]
      diag(W) <- 0
      if (!weights.scaled) W <- W/rowSums(W)
    }
  }  

  features <- features %||% rownames(object)
  features <- intersect(rownames(object),features)
  
  tab <- object[[assay]]@meta.features
  
  if (block.name %ni% colnames(tab)) {
    stop(paste0("No block.name found in the feature table of assay ", assay, ". Run LoadEPAnno first."))
  }
  
  tab <- tab[tab[[block.name]] != ".",] # skip unannotated records  
  blocks <- names(which(table(tab[[block.name]]) >= min.features.per.block))
  tab <- subset(tab, tab[[block.name]] %in% blocks)
  
  features <- intersect(rownames(tab), features)
  blocks <- unique(tab[features,][[block.name]])
  
  ## all.features <- rownames(tab)
  if (block.assay %ni% names(obj)) {
    message("No assay found in object, try to aggregate features in blocks to build a new assay.")
    if (slot != "counts") {
      warning("Automatic set slot to \"counts\", because new assay requires sum up counts by block.")
      slot <- "counts"
    }

    x <- GetAssayData(object, assay = assay, slot = slot)

    # cell sizes
    cs <- colSums(x)

    x <- x[rownames(tab),]
    x <- as(x, "TsparseMatrix")
    
    # Aggregate features in the same block
    y <- sparseMatrix(i = match(tab[[block.name]][x@i+1], blocks),
                      j = x@j+1,
                      x = x@x)
    
    rownames(y) <- blocks
    colnames(y) <- cells
    
    object[[block.assay]] <- CreateAssayObject(counts = y, assay = block.assay)

    if (alter.splice.mode) {
      ## expand matrix of block features for calculating with EP matrix
      y <- y[tab[[block.name]],]
      rownames(y) <- rownames(x)
      y <- y - x
    }
    
    x <- log1p(t(t(x)/cs) * scale.factor)
    y <- log1p(t(t(y)/cs) * scale.factor)
    
    SetAssayData(object = object, slot = "data", new.data = y, assay=block.assay)

  } else {
    if (slot != "data") {
      warning("Reset block.assay to \"data\".")
      slot <- "data"
    }

    if (alter.splice.mode) {
      warning("Notice: alter.splice.mode only can be enabled when aggregate counts.")
      alter.splice.mode <- FALSE
    }
    
    old.assay <- DefaultAssay(object)
    DefaultAssay(object) <- block.assay
    blocks <- intersect(blocks, rownames(object))
    tab <- subset(tab, tab[[block.name]] %in% blocks)    
    x <- GetAssayData(object, assay = assay, slot = slot)[rownames(tab),cells]
    y <- GetAssayData(object, assay = block.assay, slot = slot)[blocks,cells]

    DefaultAssay(object) <- old.assay
  }

  x <- x %*% W
  y <- y %*% W
  
  x <- t(scale(t(x)))
  y <- t(scale(t(y)))

  if (!alter.splice.mode) {
    ## expand matrix of block features for calculating with EP matrix
    y <- y[tab[[block.name]],]
  }
  
  # Calcuate Co-dispersion coefficient
  # r = sum(X*Y)/sqrt(sum(X*X)*sum(Y*Y))
  rs <- rowSums(x * y)
  sq <- sqrt(rowSums(x * x)*rowSums(y * y))
  
  rm(x)
  rm(y)

  gc()

  ad <- rs/sq

  tab <- object[[assay]]@meta.features
  tab[[name]] <- ad[rownames(object)]

  object[[assay]]@meta.features <- tab

  object
}

#'
#' @export
RunTwoAssayCorr <- function(object = NULL,
                            assay1 = "RNA",
                            assay2 = "EP",
                            features1 = VariableFeatures(object, assay = RNA),
                            features2 = AutoCorrFeatures(object, assay = assay2),
                            cells = NULL,
                            slot = "data",
                            weights = NULL,
                            weights.scaled = FALSE,
                            reduction = "pca",
                            dims=NULL,
                            k.nn = 20,
                            kernel.method = "dist",
                            verbose = TRUE
                            )
{
  cells <- cells %||% colnames(object)
  cells <- intersect(colnames(object), cells)

  # Make weights
  if (is.null(weights)) {
    W <- GetWeights(object = object, reduction = reduction, dims=dims, k.nn = k.nn,
                    kernel.method = kernel.method,
                    cells = cells)
  } else {
    dims <- dim(weights)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
    if (dims[1] != length(cells)) {
      cells <- intersect(colnames(weights),cells)
      W <- weights[cells, cells]
      diag(W) <- 0
      if (!weights.scaled) W <- W/rowSums(W)
    }
  }  

  if (is.null(features1)) {
    stop("Features1 is empty.")
  }

  if (is.null(features2)) {
    stop("Features2 is empty.")
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
  
  x <- x %*% W
  y <- y %*% W
  
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
aggregateGeneSets <- function(object = NULL, name = NULL, features = NULL, assay = DefaultAssay(object), scale.factor = 1e4, scaled=TRUE)
{
  if (is.null(name)) {
    stop("Set up name for the set.")
  }

  if (is.null(features)) {
    stop("Set up features for the set.")
  }
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
