anno_vcf <- function(chr, start, end, ref, alt, strand, vcf, tags, check.alt.only)
{
  sl <- .Call("anno_vcf", chr, start, end, ref, alt, strand, vcf, tags, check.alt.only)

  sl
}
#' @export
varanno <- function(chr = NULL, start = NULL, end = NULL, ref = NULL, alt = NULL, strand = NULL, vcf = NULL, tags = NULL, check.alt.only = FALSE) {
  if (is.null(chr)) stop("No chr specified.")
  if (is.null(start)) stop("No start specified.")

  df <- data.frame(chr = chr, start = start)
  if (!is.null(vcf) & !is.null(tags)) {
    sl <- anno_vcf(chr = chr, start = start, end = end, ref = ref, alt = alt, strand = strand, vcf = vcf, tags = tags, check.alt.only = check.alt.only)
    if (is.null(sl)) error("Failed to annotate VCF.")
    n <- length(tags)
    for (i in 1:n) {
      df[[tags[[i]]]] <- sl[[i]]
    }
  }
  df
}

anno_gene <- function(chr = NULL, start = NULL, end = NULL, ref = NULL, alt = NULL, strand = NULL, gtf = NULL) {
  sl <- .Call("anno_gene", chr, start, end, ref, alt, strand, gtf)
  sl
}

#' @export
annoVAR <- function(object = NULL, assay = NULL, gtf = NULL, vcf = NULL, tags = NULL, check.alt.only = FALSE)#, adjust.af = FALSE)
{
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay

  df <- object[[assay]][[]]
  
  if (length(intersect(c("chr","start","ref","alt"), colnames(df))) != 4) {
    message("Parse names ..")
    object <- ParseVAR(object)
    df <- object[[assay]][[]]
  }

  if (!is.null(tags) & !is.null(vcf)) {
    df0 <- varanno(chr=df$chr, start=as.integer(df$start), ref=df$ref, alt=df$alt, vcf = vcf, tags = tags, check.alt.only = check.alt.only)
    for (tag in tags) {
      object[[assay]][[tag]] <- df0[[tag]]
    }
  }

  if (!is.null(gtf)) {
    df0 <- anno_gene(chr = df$chr, start = as.integer(df$start), ref = df$ref, alt = df$alt, strand = df$strand, gtf = gtf)
    object[[assay]][["gene_name"]] <- df0[[1]]
    object[[assay]][["type"]] <- df0[[2]]
  }
  DefaultAssay(object) <- old.assay
  object
}
