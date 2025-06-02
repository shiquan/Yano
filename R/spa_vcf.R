## EXPERIMENTAL
vcf_sample_names <- function(vcf)
{
  .Call("vcf_sample_names", normalizePath(vcf))
}

RunSDA_vcf <- function(vcf = NULL,
                       reduction = NULL,
                       maf = NULL,
                       prune.SNN = 1/10, threads = 0, k.param=10) {
  if (is.null(vcf)) stop("No vcf file specified")
  if (is.null(reduction)) stop("No reduction data specified")

  if (is.null(maf)) stop("maf must be set.")
  if (maf < 0 || maf > 1) stop("Not valid minor allele frequence.")
  
  
  sample_names <- vcf_sample_names(vcf)
  sample_names <- intersect(sample_names, rownames(reduction))
  if (length(sample_names) == 0) {
    stop("No shared sample names found in VCF and reduction map.")
  }

  W <- GetWeights(emb = reduction[sample_names,], prune.SNN = prune.SNN, k.param=k.param)

  threads <- getCores(threads)

  df <- .Call("spa_vcf", normalizePath(vcf), W, maf, 100, threads)
  
  df
}

vcfPlot <- function(chr = NULL, pos = NULL, ref = NULL, alt = NULL, reduction = NULL, vcf = NULL, ncol=2) {
  if (is.null(chr) || is.null(pos) || is.null(reduction) || is.null(vcf)) {
    stop("chr, pos, reduction and vcf file are required.")
  }

  if (ncol(reduction) < 2) {
    stop("reduction required at least 2 dims.")
  }
  
  df <- vcf_query_value(vcf, chr, pos, ref)
  if (is.null(df)) {
    stop("No record found.")
  }
  
  alt <- alt %||% colnames(df)

  samples <- intersect(rownames(reduction), rownames(df))
  
  rd <- as.data.frame(reduction[samples,1:2])
  
  xd <- colnames(rd)

  data <- cbind(rd, df[samples,])

  if (length(alt) == 1) {
    ncol <- 1
  }

  plist <- lapply(alt, function(xx) {
    p <- ggplot(data) + geom_point(aes_string(x = xd[1], y = xd[2], color = xx))
    p 
  })

  gridExtra::grid.arrange(plist, ncol=ncol)
}

vcf_query_value <- function(vcf = NULL, chr = NULL, pos = NULL, ref = NULL) {
  if (is.null(vcf) || is.null(chr) || is.null(pos)) {
    stop("chr, pos and vcf file are required.")
  }

  df <- .Call("vcf_query_value", vcf, chr, pos, ref)
  df
}
