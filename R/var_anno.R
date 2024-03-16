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

anno_vcf <-  function(chr = NULL, start = NULL, end = NULL, ref = NULL, alt = NULL, strand = NULL, vcf = NULL, tags = NULL, adjust.af = FALSE) {
  sl <- .Call("anno_vcf", chr, start, end, ref, alt, strand, normalizePath(vcf), tags)
  sl
}

#' @export
EATanno <- function(object = NULL, assay = NULL, gtf = NULL, vcf = NULL, tags = NULL, adjust.af = FALSE)
{
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(obj)
  DefaultAssay(obj) <- assay
  if (length(intersect(c("chr","start","ref","alt"), colnames(df))) != 4) {
    message("Parse names ..")
    obj <- ParseVAR(obj)
  }

  df <- object[[assay]][[]]

  if (!is.null(tags) & !is.null(vcf)) {
    df0 <- varanno(chr=df$chr, start=df$start, ref=df$ref, alt=df$alt, vcf = vcf, tags = tags)
    for (tag in tags) {
      object[[assay]][[tag]] <- df0[[tag]]
    }

    if (isTRUE(adjust.af) & length(tags) == 1) {
      rownames(df0) <- rownames(df)
      subset(df0, is.na(df0[[tags]]) & ref == alt) %>% rownames -> sel
      df0[sel, tags] <- 1

      subset(df0, is.na(df0[[tags]]) & ref != alt) %>% rownames -> sel
      df0[sel, tags] <- 0
      object[[assay]][[tags]] <- df0[[tags]]
    }
  }
  DefaultAssay(object) <- old.assay
  object
}
