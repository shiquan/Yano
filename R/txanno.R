#' @title annoTX
#' @description Annotate gene name and genomic locations for transcript id.
#' @param object Seurat object.
#' @param assay Work assay.
#' @param gtf GTF object, load by gtf2db.
#' @param gene.name Tag name for gene name. Will be created after annotation. Default is "gene_name". 
#' @return Annotated Seurat object with annotations.
#' @export
annoTX <- function(object = NULL, assay = NULL, gtf = NULL, gene.name = "gene_name")
{
  assay <- assay %||% DefaultAssay(object)

  if (notGTF(gtf)) {
    stop("Not set a properly gtf database. Load the GTF file by gtf2db first.")
  }
  
  tx <- rownames(object[[assay]])
  sl <- .Call("tx2gene", tx, gtf)
  ug <- unique(sl[[5]])
  if (length(ug) == 1 && ug[1] == '.') {
    warnings("No genes found. Make sure you use the transcript assay and right GTF database.")
  }

  sl0 <- sl[[1]]
  fn <- rownames(object[[assay]])
  names(sl0) <- fn
  object[[assay]][["chr"]] <- sl0
  sl0 <- sl[[2]]
  names(sl0) <- fn
  object[[assay]][["start"]] <- sl0
  sl0 <- sl[[3]]
  names(sl0) <- fn
  object[[assay]][["end"]] <- sl0
  sl0 <- sl[[4]]
  names(sl0) <- fn
  object[[assay]][["strand"]] <- sl0
  sl0 <- sl[[5]]
  names(sl0) <- fn
  object[[assay]][[gene.name]] <- sl[[5]]
  
  object <- LogSeuratCommand(object)
  object
}
