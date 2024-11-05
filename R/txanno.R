#' @title tx2gene
#' @description Annotate gene name for transcript id.
#' @param object Seurat object.
#' @param assay Work assay.
#' @param gtf GTF object, load by gtf2db.
#' @param gene.name Tag name for gene name. Default is "gene_name".
#' @return Annotated Seurat object with gene_name tag
#' @export
tx2gene <- function(object = NULL, assay = NULL, gtf = NULL, gene.name = "gene_name")
{
  assay <- assay %||% DefaultAssay(object)

  if (notGTF(gtf)) {
    stop("Not set a properly gtf database. Load the GTF file by gtf2db first.")
  }
  
  tx <- rownames(object[[assay]])
  sl <- .Call("tx2gene", tx, gtf)
  ug <- unique(sl[[5]])
  if (length(ug) == 1 && ug[1] == '.') {
    warnings("No genes found. Make sure you use the transcript assay and right gtf database.")
  }
  object[[assay]][["chr"]] <- sl[[1]]
  object[[assay]][["start"]] <- sl[[2]]
  object[[assay]][["end"]] <- sl[[3]]
  object[[assay]][["strand"]] <- sl[[4]]
  object[[assay]][[gene.name]] <- sl[[5]]
  object <- LogSeuratCommand(object)
  object
}
