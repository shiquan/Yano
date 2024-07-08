anno_vcf <- function(chr, start, end, ref, alt, strand, vcf, tags, check.alt.only)
{
  sl <- .Call("anno_vcf", chr, start, end, ref, alt, strand, vcf, tags, check.alt.only)

  sl
}
#' @title varanno
#' @description Annotate genetic variants.
#' @param chr Vector of chromosome names.
#' @param start Vector of start positions.
#' @param end Vector of end positions.
#' @param ref Vector of reference alleles.
#' @param alt Vector of alternative alleles.
#' @param strand Vector of strands.
#' @param vcf VCF database. Should be indexed with `bcftools index` at first.
#' @param tags  Vector of tags to annotate. Require VCF database specified, and tags should be well formated in the VCF header.
#' @param check.alt.only Only annotate records for alternative allele (non-ref allele). Default is FASLE.
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

#' @title annoVAR
#' @description Annotate gene region and other tags for genetic variants. VCF database can be specified, and the type of VCF tags will be consider. See VCF specification and manual page for details <https://shiquan.github.io/Yano.html>.
#' @param object Seurat object.
#' @param assay Work assay.
#' @param gtf GTF object, load by gtf2db.
#' @param vcf VCF database. Should be indexed with `bcftools index` at first.
#' @param tags Vector of tags to annotate. Require VCF database specified, and tags should be well formated in the VCF header.
#' @param check.alt.only Only annotate records for alternative allele (non-ref allele). Default is FASLE.
#' @return Annotated Seurat object.
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
