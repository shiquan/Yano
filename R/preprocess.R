#' @export
setMethod(f = "QuickRecipe0",
          signature = signature(counts = "SMatrix"),
          definition = function(counts = NULL, min.cells = 20, min.features = 200,
                                assay = NULL, verbose = TRUE
                                ) {

            assay <- assay %||% "RNA"

            if (isTRUE(verbose)) {
              message(paste0("object <- CreateSeuratObject(counts = counts, min.cells = ", min.cells, ", min.features = ",min.features, ", assay = \"", assay, "\")"))
            }
            counts <- CreateSeuratObject(counts = counts, min.cells = min.cells,
                                         min.features = min.features, assay = assay)

            
            return(QuickRecipe0(counts, assay = assay, verbose = verbose))
          })

#' @export
setMethod(f = "QuickRecipe0",
          signature = signature(counts = "Seurat"),
          definition = function(counts = NULL, scale.factor = 1e4,
                                assay = NULL, verbose = TRUE
                                ) {
            assay <- assay %||% DefaultAssay(counts)
            if (verbose) {
              message(paste0("Set default assay to ", assay))
            }
            DefaultAssay(counts) <- assay

            if (verbose) {
              message(paste0("object <- NormalizeData(object, normalization.method = \"LogNormalize\", scale.factor = ", scale.factor, ")"))
            } 
            
            counts <- NormalizeData(counts, normalization.method = "LogNormalize",
                                    scale.factor = scale.factor, verbose = verbose)
            counts
          })

ProcessDimReduc <- function(object = NULL, ndim=20, resolution = 0.5, nvar= 3000, features = NULL, verbose = TRUE)
{

  if (isTRUE(verbose)) {
    message(paste0("object <- FindVariableFeatures(object, selection.method = \"vst\", nfeatures = ", nvar, ")"))
  } 

  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nvar, verbose = verbose)
  
  features <- features %||% VariableFeatures(object)
  features <- intersect(features,rownames(object))

  if (isTRUE(verbose)) {
    message("object <- ScaleData(object, features =  @features)")
  } 

  object <- ScaleData(object, features = features, verbose = verbose)

  if (isTRUE(verbose)) {
    message(paste0("object <- RunPCA(object, features = @features"))
  } 

  object <- RunPCA(object, features = features,verbose = verbose)

  if (isTRUE(verbose)) {
    message(paste0("object <- FindNeighbors(object, dims = 1:", ndim,")"))
  } 

  object <- FindNeighbors(object, dims = 1:ndim, verbose = verbose)

  if (isTRUE(verbose)) {
    message(paste0("object <- FindClusters(object, resolution = ", resolution, ")"))
  } 
  object <- FindClusters(object, resolution = resolution, verbose = verbose)

  if (isTRUE(verbose)) {
    message(paste0("object <- RunUMAP(object, dims = 1:", ndim, ")"))
  } 

  object <- RunUMAP(object, dims = 1:ndim, verbose = verbose)
  object
}

#' @title QuickRecipe
#' @description Quick clust single cell gene expression matrix with Seurat pipeline
#' @param counts raw counts matrix or Seurat object.
#' @param min.cells Only compute for features in at least this many cells
#' @param min.features Only computer for cells contained at least this many features
#' @param nvar Number of high variable features selected for PCA analysis
#' @param resolution Value of the resolution parameter pass to FindClusters, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
#' @param assay Assay name. Default is 'RNA'.
#' @param ndim Use top N PCs for clustering and UMAP. Default is 20.
#' @param verbose Print log information.
#' @return Seurat object
#' @import Seurat
#' @import Matrix
#'
#' @export
setMethod(f = "QuickRecipe",
          signature = signature(counts = "Seurat"),
          definition = function(counts = NULL, min.cells = 20, min.features = 200,
                                nvar = 3000, resolution = 0.5, assay = NULL,
                                ndim = 20, verbose = TRUE
                                ) {
            
            object <- QuickRecipe0(counts=counts, 
                                   min.cells =min.cells, min.features = min.features,
                                   assay = NULL, verbose = verbose)

            ProcessDimReduc(object, ndim=ndim, resolution=resolution, nvar=nvar, verbose = verbose)
          })

setMethod(f = "QuickRecipe",
          signature = signature(counts = "SMatrix"),
          definition = function(counts = NULL, min.cells = 20, min.features = 200,
                                nvar = 3000, resolution = 0.5, assay = NULL,
                                ndim = 20, verbose = TRUE
                                ) {
            
            object <- QuickRecipe0(counts=counts,
                                   min.cells =min.cells, min.features = min.features,
                                   assay = assay, verbose = verbose)

            ProcessDimReduc(object, ndim=ndim, resolution=resolution, nvar=nvar, verbose = verbose)
          })


"%ni%" <- Negate("%in%")

#' @importFrom parallel detectCores
#' @export
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
#' @title ParseExonName
#' @description Parse chromosome, start, end, strand, and gene name from the exon name. The exon name generated by `PISA anno` is formated like chr:start-end/[+-]/gene_name.
#' @param object Seurat object.
#' @param assay Exon assay name. Default use current actived assay.
#' @examples
#' data("glbt_small")
#' DefaultAssay(glbt_small) <- "exon"
#' # Check the meta table before parsing
#' head(glbt_small[['exon']][[]])
#' 
#' glbt_small <- ParseExonName(glbt_small)
#'
#' # Now see the meta table after parsing
#' head(glbt_small[['exon']][[]])
#' 
#' @export
ParseExonName <- function(object = NULL, assay = NULL)
{
  assay <- assay %||% DefaultAssay(object)
  message(paste0("Working on assay ", assay))
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay

  rn <- rownames(object)
  sl <- .Call("parse_exon_names", rn)

  chr <- sl[[1]]
  start <- sl[[2]]  
  end <- sl[[3]]
  gene <- sl[[4]]
  strand <- sl[[5]]
  names(chr) <- names(start) <- names(end) <- names(gene) <- rn 

  object0 <- object[[assay]]
  object0[['chr']] <- chr
  object0[['start']] <- start
  object0[['end']] <- end
  object0[['gene_name']] <- gene

  if (!is.null(strand)) {
    names(strand) <- rn
    object0[['strand']] <- strand
  }
    
  object[[assay]] <- object0
  DefaultAssay(object) <- old.assay

  object <- LogSeuratCommand(object)
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

  object <- LogSeuratCommand(object)
  object
}
#' @title ParseVarName
#' @description Parse chromosome name, start, ref allele, alternative allele and strand from EAT name. EAT generated by `PISA anno` is formated like chr:pos(ref)>(alt)/[+-].
#' @param object Seurat object
#' @param assay EAT assay. Default is current actived assay.
#' @export
ParseVarName <- function(object = NULL, assay = NULL)
{
  old.assay <- DefaultAssay(object)
  assay <- assay %||% old.assay
  DefaultAssay(object) <- assay

  rn <- rownames(object)
  sl <- .Call("parse_var_names", rn)

  if (is.null(sl)) {
    message("Failed to parse variant name.")
    return(NULL)
  }
  chr <- sl[[1]]
  start <- sl[[2]]
  ref <- sl[[3]]
  alt <- sl[[4]]
  strand <- sl[[5]]
  names(chr) <- names(start) <- names(ref) <- names(alt) <- rn 

  object0 <- object[[assay]]
  object0[['chr']] <- chr
  object0[['start']] <- start
  object0[['ref']] <- ref
  object0[['alt']] <- alt


  if (is.null(strand)) {
    locs <- paste0(sl[[1]],":",sl[[2]])
  } else {
    names(strand) <- rn
    object0[['strand']] <- strand
    locs <- paste0(sl[[1]],":",sl[[2]],"/",strand)
  }
  
  names(locs) <- rn
  
  object0[['locus']] <- locs
  
  object[[assay]] <- object0
  DefaultAssay(object) <- old.assay

  object <- LogSeuratCommand(object)
  object
}
#' @title ParseBED
#' @description Parse chromosome, start, end, strand from the BED/EPT name. The BED name is generated by `PISA annotate`.
#' @param object Seurat object
#' @param assay EPT assay. Default is current actived assay.
#' @export
ParseBED <- function(object = NULL, assay = NULL)
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

  object <- LogSeuratCommand(object)
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
