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

            #RNA.assay <- intersect("RNA", names(counts))
            #RNA.assay <- RNA.assay %||% assay
            #message(paste0("Use assay ", RNA.assay, " to evaluate cell size.." ))
            #x <- GetAssayData(counts, assay = RNA.assay, slot = "counts")
            #cs <- colSums(x)
            #cells <- intersect(colnames(counts),names(which(cs > 0)))
            #cs <- cs[cells]
            #counts <- counts[,cells]
            #x <- GetAssayData(counts, assay=assay, slot = "counts")
            #x <- log1p(t(t(x)/cs) * scale.factor)
            #rownames(x) <- rownames(counts)
            #colnames(x) <- cells
            #SetAssayData(object = counts, slot = "data", new.data = x, assay= assay)
            #rm(x)
            #gc()
            counts
          })

#' @export
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

#'
#' @export
GetWeights <- function(object= NULL,                       
                       reduction = "pca",
                       dims = NULL,
                       k.nn = 5,
                       spatial = FALSE,
                       kernel.method = "average",
                       self.weight = 0,
                       scale = FALSE,
                       cells = NULL)
{
  cells <- cells %||% colnames(object)
  
  if (isTRUE(spatial)) {
    message("Build weights on tissue coordiantes")
    emb <- GetTissueCoordinates(object)
    kernel.method <- "average"
  } else {
    message(paste0("Build weights on ", reduction))
    emb <- Embeddings(object,reduction = reduction)
    if (!is.null(dims)) emb <- emb[,dims]
  }

  emb <- emb[cells,]
  
  knn.rlt <- nabor::knn(data=emb, query = emb, k=k.nn)
  if (kernel.method == "average") x = 1
  else x = 1/c(knn.rlt$nn.dists)

  ncell <- length(cells)
  W <- sparseMatrix(i = rep(c(1:ncell), k.nn), 
                    j = c(knn.rlt$nn.idx),
                    x = x,
                    dims = c(ncell, ncell))

  diag(W) <- self.weight
  
  if (scale) W <- W/rowSums(W)
  W[is.na(W)] <- 0
  colnames(W) <- cells
  rownames(W) <- cells
  W
}

#'
#' @import sparseMatrixStats
#' @import Matrix
#' @export
RunAutoCorr <- function(object = NULL, 
                        assay = DefaultAssay(object = object),
                        slot = "data",
                        spatial = FALSE,
                        scaled = FALSE,
                        weights = NULL,                              
                        scale.weight = FALSE,
                        reduction = "pca",
                        dims = NULL,
                        k.nn = 5,
                        kernel.method = "dist",
                        cells = NULL,
                        features = NULL,
                        perm = 10000,threads=4,
                        verbose = TRUE)
{
  message(paste0("Working on assay : ", assay))
  cells <- cells %||% colnames(object)
  cells <- intersect(colnames(object), cells)

  features <- features %||% rownames(object)
  features <- intersect(rownames(object),features)
  
  if (is.null(weights)) {
    W <- GetWeights(object=object, reduction=reduction,
                    dims=dims,
                    k.nn=k.nn,
                    kernel.method=kernel.method,
                    cells=cells,
                    self.weight = 0,
                    spatial=spatial)
  } else {
    dims <- dim(weights)
    if (dims[1] != dims[2]) stop("Weight matrix should be a squared matrix.")
    if (dims[1] != length(cells)) {
      cells <- intersect(colnames(weights),cells)
      W <- weights[cells, cells]
      diag(W) <- 0
      if (scale.weight) W <- W/rowSums(W)
    }
  }
  
  x0 <- GetAssayData(object, assay = assay, slot = slot)[,cells]
  features <- intersect(rownames(x0), features)
  x0 <- x0[features,]
  x0 <- as(x0, "CsparseMatrix")
  W <- as(W, "CsparseMatrix")
  
  message(paste0("Run autocorrelation test for ", length(features), " features."))
  moransi.vals <- .Call("autocorrelation_test", x0, W, TRUE, threads);

  #message("Run permutation test.")
  #moransi.vals2 <-  .Call("moransi_mc_test", x0, W, TRUE, perm, threads);

  ## message(paste0("Run Geary's C for ", length(features), " features."))
  ## gearysc.vals <- .Call("GearysC_test", x0, W);
  Ivals <- moransi.vals[[1]]
  Cvals <- moransi.vals[[2]]
  IZvals <- moransi.vals[[3]]  
  CZvals <- moransi.vals[[4]]
  #Ivals2 <- moransi.vals2[[1]]
  #Pvals <- moransi.vals2[[2]]
  names(Ivals) <- features
  names(Cvals) <- features
  names(IZvals) <- features
  names(CZvals) <- features
  #names(Ivals2) <- features
  #names(Pvals) <- features
  ## names(gearysc.vals) <- features

  object[[assay]]@meta.features[['MoransI']] <- Ivals[rownames(object)]
  object[[assay]]@meta.features[['MoransI.Z']] <- IZvals[rownames(object)]
  object[[assay]]@meta.features[['GearyC']] <- Cvals[rownames(object)]
  object[[assay]]@meta.features[['GearyC.Z']] <- CZvals[rownames(object)]
  #object[[assay]]@meta.features[['MoransI.2']] <- Ivals2[rownames(object)]
  #object[[assay]]@meta.features[['MoransI.pval']] <- Pvals[rownames(object)]
  rm(W)
  rm(moransi.vals)
  #rm(moransi.vals2)
  ## rm(gearysc.vals)
  gc()
  
  object
}

#' Set autocorrection by rank and/or score
#' @export
SetAutoCorrFeatures <- function(object = NULL,
                                moransi.min = 0,
                                top.n = 500,
                                assay = DefaultAssay(object),
                                sd = 3,
                                degree=1,
                                plot = TRUE,
                                return.plot = FALSE
                                )
{
  tab <- object[[assay]]@meta.features

  cn <- colnames(tab)
  if ("MoransI.value" %ni% cn | "MoransI.rank" %ni% cn | "coverage" %ni% cn) {
    stop("No Morans'I value found, use RunAutoCorr first.")
  }

  idx <- which(tab[["MoransI.value"]] > moransi.min & tab[["MoransI.rank"]] <= top.n)

  n1 <- length(idx)
  
  los <- loess.smooth(tab$coverage, tab$MoransI.value, degree=degree)
  func <- approxfun(los$x, los$y * sd)
  
  idx <- intersect(which(tab[['MoransI.value']] > func(tab[['coverage']])),idx)
  n2 <- length(idx)
  
  tab[["AutoCorrFeature"]] <- FALSE
  tab[idx,][["AutoCorrFeature"]] <- TRUE
  object[[assay]]@meta.features <- tab

  if (plot) {
    v <- min(tab[idx,][["MoransI.value"]])
    
    p1 <- ggplot(tab, aes(x=MoransI.value)) + geom_density(size=1) + theme_minimal()
    p1 <- p1 + geom_vline(aes(xintercept=v),color="red", linetype="dashed", size=1)
    p1 <- p1 + ggtitle(paste0("n = ",n1))
    p2 <- ggplot() + geom_point(data=tab,aes(x=coverage,y=MoransI.value),color="grey")
    p2 <- p2 + geom_line(data=as.data.frame(los),aes(x,y),size=1, color="black") + theme_minimal()
    p2 <- p2 + geom_line(data=as.data.frame(los),aes(x,y*sd),size=1, linetype="dashed", color="blue")    
    p2 <- p2 + geom_point(data=tab[idx,], aes(coverage, MoransI.value), shape=21,size=3)
    p2 <- p2 + geom_hline(aes(yintercept=v),color="red", linetype="dashed", size=1)
    p2 <- p2 + theme_minimal()
    p2 <- p2 + ggtitle(paste0("n = ",n2))
    p <- cowplot::plot_grid(p1,p2)
    print(p)
    if (isTRUE(return.plot)) {
      return(p)
    }
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
                      features = NULL,
                      cells = NULL,
                      assay = DefaultAssay(object),             
                      slot = "data",
                      scaled = FALSE,
                      weights = NULL,
                      weights.scaled = FALSE,
                      reduction = "pca",
                      dims=NULL,
                      k.nn = 5,
                      kernel.method = "average",
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

  if (is.null(weights)) {
    W <- GetWeights(object = object, reduction = reduction,
                    dims = dims,
                    k.nn = k.nn,
                    kernel.method = kernel.method,
                    self.weight = self.weight,
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

  # mtx <- mtx %*% W
  mtx <- Matrix::tcrossprod(mtx,W)

  # smoothed by weights for logcounts, then scaled 
  if (!scaled) {
    ## todo, performance improvement
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
#'
#' @importFrom data.table fread
#' @export
LoadBEDanno <- function(file = NULL, object = NULL, assay = NULL, stranded = TRUE)
{
  bed <- fread(file)[,c(1:9)]
  colnames(bed) <- c("chr","start","end","name","score","strand","n_gene","gene_name","type")
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

  assay <- assay %||% DefaultAssay(object)
  
  object[[assay]]@meta.features[['chr']] <- bed[['chr']]
  object[[assay]]@meta.features[['start']] <- bed[['start']]
  object[[assay]]@meta.features[['end']] <- bed[['end']]
  object[[assay]]@meta.features[['name']] <- bed[['name']]
  object[[assay]]@meta.features[['strand']] <- bed[['strand']]
  object[[assay]]@meta.features[['n_gene']] <- bed[['n_gene']]
  object[[assay]]@meta.features[['gene_name']] <- bed[['gene_name']]
  object[[assay]]@meta.features[['type']] <- bed[['type']]
  
  object
}

#'
#' @export
RunBlockCorr <- function(object = NULL,
                         block.name = "gene_name",
                         assay = NULL,
                         name = NULL,
                         features = NULL,
                         sensitive.mode = FALSE,
                         block.assay = NULL,
                         block.assay.replace = FALSE,
                         cells = NULL,
                         feature.types = c("exon","exonintron","intron","multiexons","utr3","utr5"),
                         min.features.per.block = 2,
                         scale.factor = 1e4,
                         weights = NULL,
                         weights.scaled = FALSE,
                         reduction = "pca",
                         dims = NULL,
                         k.nn = 5,
                         kernel.method = "average",
                         self.weight = 1,
                         #keep.matrix = FALSE,
                         perm.test = FALSE,
                         perm=1000,
                         threads = 1,
                         verbose = TRUE
                         )
{
  cells <- cells %||% colnames(object)
  cells <- intersect(colnames(object), cells)

  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))
  
  name <- name %||% block.name

  features <- features %||% AutoCorrFeatures(object)
  features <- intersect(rownames(object),features)

  message(paste0("Working on ", length(features), " features."))
  
  # Make weights
  if (is.null(weights)) {
    W <- GetWeights(object = object,
                    reduction = reduction,
                    dims=dims,
                    k.nn = k.nn,
                    kernel.method = kernel.method,
                    cells = cells,
                    self.weight = self.weight,
                    scale=TRUE)
                   
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

  tab <- object[[assay]]@meta.features

  if (block.name %ni% colnames(tab)) {
    stop(paste0("No block.name found in the feature table of assay ", assay, ". Run LoadEPTanno first."))
  }
  
  tab <- tab[tab[[block.name]] != ".",] # skip unannotated records
  if ("type" %in% colnames(tab)) {
    tab <- subset(tab, type %in% feature.types)
  }
  
  blocks <- names(which(table(tab[[block.name]]) >= min.features.per.block))
  
  tab0 <- tab[features,]
  blocks <- intersect(unique(tab0[[block.name]]), blocks)

  message(paste0("Processing ", length(blocks), " blocks.."))
  tab <- subset(tab, tab[[block.name]] %in% blocks)
  
  if (length(features) == 0) {
    stop("No features found.")
  }

  block.assay <- block.assay %||% "tmp.assay"
  ## blocks <- unique(blocks)

  x <- GetAssayData(object, assay = assay, slot = "counts")
  
  # cell sizes
  cs <- colSums(x)
  cells <- intersect(cells,names(which(cs > 0)))
  W <- W[cells, cells]
  cs <- cs[cells]
  ncell <- length(cells)

  ## all.features <- rownames(tab)
  if (block.assay %ni% names(obj) || block.assay.replace) {
    message("Aggregate counts..")
    #if (slot != "counts") {
     # warning("Automatic set slot to \"counts\", because new assay requires sum up counts by block.")
      #slot <- "counts"
    #}
  
    ## if (ncell != ncol(object) && keep.matrix) {
    ##   warnings("Inconsistance cells, set keep.matrix to FALSE.")
    ##   keep.matrix <- FALSE
    ## }

    #x <- GetAssayData(object, assay = assay, slot = "counts")
    x <- x[rownames(tab), cells]
    x <- as(x, "TsparseMatrix")
    
    # Aggregate features in the same block
    y <- sparseMatrix(i = match(tab[[block.name]][x@i+1], blocks),
                      j = x@j+1,
                      x = x@x, dims=c(length(blocks), length(cells)))
    
    rownames(y) <- blocks
    colnames(y) <- cells

    ## if (keep.matrix) {
    ##   object[[block.assay]] <- CreateAssayObject(counts = y, assay = block.assay)
    ## }

    y <- y[tab[[block.name]],]
    rownames(y) <- rownames(x)
    
    if (sensitive.mode) {
      ## expand matrix of block features for calculating with EPT matrix
      y <- y - x
    }

  } else {

    message(paste0("Trying to retrieve data from assay ", block.assay,".."))
    
    if (isTRUE(sensitive.mode)) {
      # message("Notice: sensitive.mode only can be enabled when aggregate counts.")
      sensitive.mode <- FALSE
    }
    
    old.assay <- DefaultAssay(object)
    DefaultAssay(object) <- block.assay
    blocks <- intersect(blocks, rownames(object))
    tab <- subset(tab, tab[[block.name]] %in% blocks)

    x <- x[rownames(tab), cells]
    #x <- GetAssayData(object, assay = assay, slot = "counts")[rownames(tab),cells]
    y <- GetAssayData(object, assay = block.assay, slot = "counts")[blocks,cells]    
    DefaultAssay(object) <- old.assay

    y <- y[tab[[block.name]],]
    rownames(y) <- rownames(tab)

  }

  feature.names <- rownames(x)
  fc <- log1p(rowMeans(x)) - log1p(rowMeans(y))
  names(fc) <- feature.names

  x <- log1p(t(t(x)/cs) * scale.factor)
  y <- log1p(t(t(y)/cs) * scale.factor)
  
  #if (keep.matrix) {
  #  SetAssayData(object = object, slot = "data", new.data = y, assay=block.assay)
  #}
  gc()
  message("Smooth data..")
  rownames(y) <- feature.names

  if (isTRUE(perm.test)) {
    ta <- .Call("E_test", x, y, W, perm, threads);
  } else {
    ta <- .Call("E_test", x, y, W, 0, threads);
  }
 
  gc()

  Lx <- ta[[1]]
  Ly <- ta[[2]]
  r <- ta[[3]]
  e <- ta[[4]]

  names(Lx) <- feature.names
  names(Ly) <- feature.names
  names(r) <- feature.names
  names(e) <- feature.names

  tab <- object[[assay]]@meta.features
  tab[[paste0(name, ".e.coef")]] <- e[rownames(object)]
  tab[[paste0(name, ".r")]] <- r[rownames(object)]
  tab[[paste0(name, ".Lx")]] <- Lx[rownames(object)]
  tab[[paste0(name, ".Ly")]] <- Ly[rownames(object)]
  tab[[paste0(name, ".fc")]] <- fc[rownames(object)]

  if (perm.test) {
    pval <- ta[[5]]
    names(pval) <- feature.names
    tab[[paste0(name, ".pval")]] <- pval[rownames(object)]
  }
  
  object[[assay]]@meta.features <- tab

  rm(ta)
  gc()
  
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
                            weights = NULL,
                            weights.scaled = FALSE,
                            reduction = "pca",
                            dims=NULL,
                            k.nn = 5,
                            kernel.method = "average",
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
  if (is.null(weights)) {
    W <- GetWeights(object = object, reduction = reduction, dims=dims, k.nn = k.nn,
                    kernel.method = kernel.method, self.weight = self.weight,
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
  
  #x <- x %*% W
  #y <- y %*% W

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
