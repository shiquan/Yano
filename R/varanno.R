anno_vcf <- function(chr, start, end, ref, alt, strand, vcf, tags, check.alt.only)
{
  sl <- .Call("anno_vcf", chr, start, end, ref, alt, strand, vcf, tags, check.alt.only)

  sl
}
#' @title varanno
#' @description Annotate genetic variants with VCF databases.
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

#' @export
anno_gene <- function(chr = NULL, start = NULL, end = NULL, ref = NULL, alt = NULL, strand = NULL, gtf = NULL) {
  sl <- .Call("anno_gene", chr, start, end, ref, alt, strand, gtf)
  sl
}
#' @export
anno_conseq <- function(chr = NULL, pos = NULL, ref = NULL, alt = NULL, strand = NULL, gtf = NULL, fasta=NULL, debug=FALSE) {
  sl <- .Call("anno_conseq", chr, as.integer(pos), ref, alt, strand, gtf, fasta, debug)
  sl
}

#' @title annoVAR
#' @description Annotate genetic variants with preload GTF and/or VCF databases. Will generate gene region and other tags.
#' @param object Seurat object.
#' @param assay Work assay.
#' @param gtf GTF object, load by gtf2db.
#' @param vcf VCF database. Should be indexed with `bcftools index` at first.
#' @param tags Vector of tags to annotate. Require VCF database specified, and tags should be well formated in the VCF header. See VCF specification (https://samtools.github.io/hts-specs/VCFv4.2.pdf) for details.
#' @param check.alt.only Only annotate records for alternative allele (non-ref allele). Default is FASLE.
#' @param chr Colname name for chromsome, default is "chr".
#' @param fasta Genome reference in FATSA format, should be indexed with `samtools faidx` first.
#' @return Annotated Seurat object.
#' @export
annoVAR <- function(object = NULL, assay = NULL, gtf = NULL, vcf = NULL, tags = NULL, check.alt.only = FALSE, chr = "chr", fasta = NULL)
{
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay

  df <- object[[assay]][[]]

  if (length(intersect(c("chr","start","ref","alt"), colnames(df))) != 3) {
    message("Parse names ..")
    object <- ParseVarName(object)
    df <- object[[assay]][[]]
  }

  if (chr %ni% colnames(df)) {
    stop("No chr found at feature meta table.")
  }
  if (!is.null(tags) & !is.null(vcf)) {
    df0 <- varanno(chr=df[[chr]], start=as.integer(df$start), ref=df$ref, alt=df$alt, vcf = vcf, tags = tags, check.alt.only = check.alt.only)
    for (tag in tags) {
      object[[assay]][[tag]] <- df0[[tag]]
    }
  }

  if (!is.null(gtf)) {
    if (is.null(fasta)) {
      df0 <- anno_gene(chr = df[[chr]], start = as.integer(df$start), ref = df$ref, alt = df$alt, strand = df$strand, gtf = gtf)
      object[[assay]][["gene_name"]] <- df0[[1]]
      object[[assay]][["type"]] <- df0[[2]]
    } else {
      df0 <- anno_conseq(chr = df[[chr]], pos = as.integer(df$start), ref = df$ref, alt = df$alt, strand = df$strand, gtf = gtf, fasta=fasta)
      if (length(df0) == 2) {
        object[[assay]][["gene_name"]] <- df0[[1]]
        object[[assay]][["consequence"]] <- df0[[2]]
      } else {
        stop(df0)
      }
    }
  }
  DefaultAssay(object) <- old.assay

  object <- LogSeuratCommand(object)
  object
}
