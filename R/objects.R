#'
#' @rdname CreateEptAssay
#' @import data.table
#'
#' @export
#' 
CreateEptAssay <- function(counts,
                           assay = 'EPT',
                           meta.data = NULL,
                           min.cells = 10,
                           row.names = NULL,
                           sep = c("-","-","-"),
                           annotated.bed)
{
  if (missing(annotated.bed)) stop("No annotated bed file set. Run `PISA annobed` first.")
  ann <- read.table(annotated.bed)
  colnames(ann) <- c("chr","start","end","name","score","strand","n_gene","gene_name","type","nearest_gene","distance")
  ann$name <- paste(ann$chr, ann$start, ann$end, ann$strand, sep="-")  
  rownames(ann) <- ann$name
  rownames(counts) <- gsub("_","-",rownames(counts))
  nm <- intersect(rownames(counts),ann$name)
  ann <- ann[nm,]
  if (length(nm)==0) stop("Empty EPTs.")
  
  obj <- CreateAssayObject(
    counts = counts[nm,],
    assay = assay,
    min.cells = min.cells)

  obj@meta.features = ann[rownames(obj),]
  
  return(obj)
}
