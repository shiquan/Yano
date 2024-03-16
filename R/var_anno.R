#' @export
varanno <- function(chr = NULL, start = NULL, end = NULL, ref = NULL, alt = NULL, strand = NULL, gtf = NULL, vcf = NULL, tags = NULL) {
  if (is.null(chr)) stop("No chr specified.")
  if (is.null(start)) stop("No start specified.")

  df <- data.frame(chr = chr, start = start)
  if (!is.null(end)) df[['end']] <- end
  if (!is.null(ref)) df[['ref']] <- ref
  if (!is.null(alt)) df[['alt']] <- alt
  if (!is.null(strand)) df[['strand']] <- strand
  
  if (!is.null(gtf)) {
    sl <- anno_hgvs(chr, start, end, ref, alt, strand, gtf)
    if (is.null(sl)) error("Failed to annotate genes.")
    
    df[['gene']] <- sl[[1]]
    df[['type']] <- sl[[2]]
  }

  if (!is.null(vcf) & !is.null(tags)) {
    sl <- anno_vcf(chr, start, end, ref, alt, strand, vcf, tags)
    if (is.null(sl)) error("Failed to annotate VCF.")
    n <- length(tags)
    for (i in 1:n) {
      df[[tags[[i]]]] <- sl[[i]]
    }
  }
  df
}

anno_hgvs <- function(chr = NULL, start = NULL, end = NULL, ref = NULL, alt = NULL, strand = NULL, gtf = NULL) {
  sl <- .Call("anno_hgvs", chr, start, end, ref, alt, strand, normalizePath(gtf))
  sl
}

anno_vcf <-  function(chr = NULL, start = NULL, end = NULL, ref = NULL, alt = NULL, strand = NULL, vcf = NULL, tags = NULL) {
  sl <- .Call("anno_vcf", chr, start, end, ref, alt, strand, normalizePath(vcf), tags)
  sl
}

EATanno <- function(object = NULL, assay = NULL, gtf = NULL, vcf = NULL, tags = NULL)
{
  assay <- assay %||% DefaultAssay(object)
  df <- object[[assay]][[]]

  nm <- rownames(df)
  
  
}
