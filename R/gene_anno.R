#' @title annoGene
#' @description Annotate genomic locations for gene names.
#' @param object Seurat object.
#' @param assay Work assay.
#' @param gtf GTF object, load by gtf2db.
#' @param gene.name Tag name for gene name. If not set, use rownames instead and will create gene_name tag in the meta table.
#' @return Annotated Seurat object with annotations.
#' @export
annoGene <- function(object = NULL, assay = NULL, gtf = NULL, gene.name = NULL)
{
  assay <- assay %||% DefaultAssay(object)

  if (notGTF(gtf)) {
    stop("Not set a properly gtf database. Load the GTF file by gtf2db first.")
  }

  genes <- rownames(object[[assay]])
  if (!is.null(gene.name)) {
    genes <- object[[assay]][[gene.name]]
  } else {
    names(genes) <- genes
    object[[assay]][["gene_name"]] <- genes
  }
  
  sl <- .Call("gene_locs", genes, gtf)
  
  if (class(sl) == "character") {
    warnings("No genes found. Make sure you use the gene assay and right GTF database.")
    return(object)
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

  object <- LogSeuratCommand(object)
  object
}
