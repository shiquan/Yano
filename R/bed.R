bedanno <- function(chr = NULL, start = NULL, end = NULL, strand = NULL, gtf = NULL, upstream = 1000, downstream = 1000)
{
  sl <- .Call("anno_bed", chr, start, end, strand, gtf, upstream, downstream)
  sl
}

#' @export
annoBED <- function(object = NULL, assay = NULL, gtf = NULL)
{
  if (notGTF(gtf)) stop("GTF is not specific, use gtf2db load GTF file first.")
  
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay

  df <- object[[assay]][[]]
  
  if (length(intersect(c("chr","start","end","strand"), colnames(df))) != 4) {
    message("Parse names ..")
    object <- ParseBED(object)
    df <- object[[assay]][[]]
  }

  df0 <- bedanno(chr=df$chr, start=as.integer(df$start), end=df$end, strand = df$strand, gtf=gtf)

  object[[assay]][["gene_name"]] <- df0[[1]]
  object[[assay]][["type"]] <- df0[[2]]

  DefaultAssay(object) <- old.assay
  object
}
