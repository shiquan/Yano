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

#'@importFrom parallel detectCores
getCores <- function(threads = 0)
{
  check <- .Call("openmp_support")
  if (isTRUE(check)) {
    if (threads > 0) return(threads)
    threads <- detectCores() - 1
    if (threads > 1) return(threads)
    return(1)
  }
  if (threads > 1) {
    message("No openmp support, multithreads disabled.")
  }
  return(1)
}

ParseExonName.Assay <- function(object)
{
  nm <- rownames(object)

  chr <- gsub("(.*):(.*)-(.*)/(.)/(.*)","\\1",nm)
  start <- as.integer(gsub("(.*):(.*)-(.*)/(.)/(.*)","\\2",nm))
  end <- as.integer(gsub("(.*):(.*)-(.*)/(.)/(.*)","\\3",nm))
  strand <- gsub("(.*):(.*)-(.*)/(.)/(.*)","\\4",nm)
  gene_name <- gsub("(.*):(.*)-(.*)/(.)/(.*)","\\5",nm)
  names(chr) <- nm
  names(start) <- nm
  names(end) <- nm
  names(strand) <- nm
  names(gene_name) <- nm
  
  object[['chr']] <- chr
  object[['start']] <- start
  object[['end']] <- end
  object[['strand']] <- strand 
  object[['gene_name']] <- gene_name
  
  return(object)
}
#' @export
ParseExonName <- function(object = NULL, assay = NULL)
{
  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  object[[assay]] <- ParseExonName.Assay(object[[assay]])
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
  object0 <- object[[assay]]
  object0[[colnames(bed)]] <- bed
  object[[assay]] <- object0

  DefaultAssay(object) <- old.assay

  object
}

#'@export
ParseVAR <- function(object = NULL, assay = NULL, ignore.strand = FALSE)
{
  old.assay <- DefaultAssay(object)
  assay <- assay %||% old.assay
  DefaultAssay(object) <- assay
  rn <- rownames(object)

  sl <- .Call("parse_var_names", rn)

  if (is.null(sl)) return(NULL);
  chr <- sl[[1]]
  start <- sl[[2]]
  ref <- sl[[3]]
  alt <- sl[[4]]
  strand <- sl[[5]]
  names(chr) <- names(start) <- names(ref) <- names(alt) <- names(strand) <- rn 

  object0 <- object[[assay]]
  object0[['chr']] <- chr
  object0[['start']] <- start
  object0[['ref']] <- ref
  object0[['alt']] <- alt
  object0[['strand']] <- strand

  locs <- paste0(sl[[1]],":",sl[[2]],"/",sl[[5]])
  names(locs) <- rn
  
  object0[['locus']] <- locs
  
  object[[assay]] <- object0
  DefaultAssay(object) <- old.assay

  object
}

#'@export
ParseBED <- function(object = NULL, assay = NULL, ignore.strand = FALSE)
{
  old.assay <- DefaultAssay(object)
  assay <- assay %||% old.assay
  DefaultAssay(object) <- assay
  rn <- rownames(object)

  sl <- .Call("parse_bed_names", rn)

  if (is.null(sl)) return(NULL);
  chr <- sl[[1]]
  start <- sl[[2]]
  end <- sl[[3]]
  strand <- sl[[4]]

  names(chr) <- names(start) <- names(end) <- names(strand) <- rn 

  object0 <- object[[assay]]
  object0[['chr']] <- chr
  object0[['start']] <- start
  object0[['end']] <- end
  object0[['strand']] <- strand
  
  object[[assay]] <- object0
  DefaultAssay(object) <- old.assay

  object
}

#' @export
RenameVARs <- function(counts = NULL, strategy = 1) {
  rn <- rownames(counts)
  check <- grep("(.*:[0-9]+)([ACGT=>]*).*", rn, invert = TRUE)
  
  if (length(check) > 3) {
    stop(paste0("Unrecongised names :", paste(rn[check[1:3]], sep=",")), " ...")
  } else if (length(check) > 0) {
    stop(paste0("Unrecongised names :", paste(rn[check[1:3]], sep=",")))
  }
  
  if (strategy == 2) {
    nm <- unlist(lapply(rn, function(x) {
      s <- length(grep("/[+-]$", x))
      str <- length(grep("=", x))
      if (str == 0) {
        if (s >0) {
          gsub("(.*:[0-9]+[ACGT]*>).*/([-+])", "\\1B/\\2", x)
        } else {
          gsub("(.*:[0-9]+[ACGT]*>).*", "\\1B", x)
        }
      } else {
        x
      }
    }))

    nm0 <- unique(nm)

    idx <- match(nm, nm0)
    #counts <- as(counts, "dgCMatrix")
    x <- Matrix::sparseMatrix(i = idx[counts@i+1], p = counts@p, x= counts@x)
    rownames(x) <- nm0
    colnames(x) <- colnames(counts)
    return(x)
  } else if (strategy == 3) {

    locs <- unlist(lapply(rn, function(x) {
      s <- length(grep("/[+-]$", x))
      if (s > 0) {
        gsub("(.*:[0-9]+)([ACGT=>]*).*/([-+])", "\\1/\\3", x)
      } else {
        gsub("(.*:[0-9]+)([ACGT=>]*).*", "\\1", x)
      }
    }))
    
    idx <- names(which(table(locs) > 1))    
    sel <- rn[which(locs %in% idx)]
    
    counts <- counts[sel,]

    rm <- rowMeans(counts)
    #names(locs) <- rn
    nm <- unlist(lapply (idx, function(x) {
      sel0 <- rn[which(locs %in% x)]
      s <- length(grep("/[+-]$", x))
      if (s > 0) {
        alleleA <- gsub("(.*)(/.*)","\\1A\\2",x)
        alleleB <- gsub("(.*)(/.*)","\\1B\\2",x)
      } else {
        alleleA <- paste0(x,"A")
        alleleB <- paste0(x,"B")
      }
      new <- rep(alleleB, length(sel0))

      mi <- which.max(rm[sel0])[1]
      new[mi] <- alleleA
      names(new) <- sel0
      new
    }))
    
    nm0 <- unique(nm)
    idx <- match(nm, nm0)
    counts <- counts[names(nm),]
    x <- Matrix::sparseMatrix(i = idx[counts@i+1], p = counts@p, x= counts@x)
    rownames(x) <- nm0
    colnames(x) <- colnames(counts)
    return(x)
    
  } else {
    stop("Only support strategy 2 or 3.")
  }
  NULL
}
