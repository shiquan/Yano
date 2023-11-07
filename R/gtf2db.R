#'@useDynLib Yano C_gtf2db
#'@export
gtf2db <- function(filename = NULL) {
  
  if (is.null(filename)) stop("No gtf file.")
  tab <- .Call("C_gtf2db", normalizePath(filename))
  gene <- tab[[1]]
  transcript <- tab[[2]]
  exon <- tab[[3]]

  anno <- SimpleList()

  anno$gene <- GRanges(seqnames = gene$chr, strand = gene$strand,
                       ranges = IRanges(start = gene$start, end=gene$end),
                       gene_name = gene$genename, gene_id = gene$geneid,
                       biotype = gene$biotype)

  anno$transcript <- GRanges(seqnames = transcript$chr, strand = transcript$strand,
                       ranges = IRanges(start = transcript$start, end=transcript$end),
                       gene_name = transcript$genename, gene_id = transcript$geneid,
                       transcript_id = transcript$txid,
                       biotype = transcript$biotype)

  anno$exon <- GRanges(seqnames = exon$chr, strand = exon$strand,
                       ranges = IRanges(start = exon$start, end=exon$end),
                       gene_name = exon$genename, gene_id = exon$geneid,
                       transcript_id = exon$txid)

  rm(tab)

  return(anno)
}
