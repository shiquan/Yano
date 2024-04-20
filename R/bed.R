bedanno <- function(chr = NULL, start = NULL, end = NULL, strand = NULL, gtf = NULL, promoter = FALSE, upstream = 1000, downstream = 0, at_upstream = 1000, at_downstream = 1000)
{
  sl <- .Call("anno_bed", chr, start, end, strand, gtf, upstream, downstream, at_upstream, at_downstream)
  sl
}

#' @export
annoBED <- function(object = NULL, assay = NULL, gtf = NULL)
{
  if (notGTF(gtf)) stop("GTF is not specific, use gtf2db load GTF file first.")
  
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(obj)
  DefaultAssay(obj) <- assay

  df <- object[[assay]][[]]
  
  if (length(intersect(c("chr","start","end","strand"), colnames(df))) != 4) {
    message("Parse names ..")
    obj <- ParseBED(obj)
  }

  df0 <- bedanno(chr=df$chr, start=as.integer(df$start), end=df$end, strand = df$strand, gtf=gtf)

  object[[assay]][["gene"]] <- df0[[1]]
  object[[assay]][["type"]] <- df0[[2]]

  DefaultAssay(object) <- old.assay
  object
}
