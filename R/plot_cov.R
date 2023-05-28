#'@useDynLib Yano depth2matrix
#'@export
bamcov <- function(bamfile = NULL, chr = NULL, start = -1, end = -1, strand = "both", split.bc = FALSE, cell.group = NULL, bin=1000, cell.tag = "CB", umi.tag = "UB")
{
  if (is.null(bamfile)) stop("No BAM.")
  if (is.null(chr)) stop("No chromosome name")
  if (start == -1) stop("Start position.")
  if (end == -1) stop("End position")

  len <- end - start
  if (len > 10000000) { ## 10M
    stop("Too large region, try to shrink it before querying.")
  }

  win <- 1
  if (len > bin*2) {
    win <- as.integer(len/bin)
  }
  
  #bamfile <- normalizePath(bamfile)
  print(bamfile)

  strand.flag <- -1
  if (strand == "forward") strand.flag <- 0
  else if(strand == "reverse") strand.flag <- 1
  else if (strand == "ignore") strand.flag <- -2

  if (!is.null(cell.group)) {
    #if (length(cell.group)) stop("cell.group should be a vector of cell groups with cell names.")    
    cell.names <- names(cell.group)
    if (is.null(cell.names)) stop("Names of cell.group should not be empty.")
    cell.group <- as.character(cell.group)
    groups <- unique(cell.group)
    group.ids <- match(cell.group, groups)
    n <- length(cell.group)
    split.bc <- TRUE    
    dlst<- .Call("depth2matrix", bamfile, chr, start, end, strand.flag, split.bc, cell.tag, umi.tag, 20, cell.names, n, group.ids, groups)
  } else {
    dlst <- .Call("depth2matrix", bamfile, chr, start, end, strand.flag, split.bc, cell.tag, umi.tag, 20, NULL, 0, NULL, NULL)
  }
  
  idx.0 <- which(dlst[[2]] == 0)
  idx.1 <- which(dlst[[2]] == 1)

  dlst[[1]] <- as.integer(dlst[[1]]/win)*win

  x <- seq(as.integer(start/win)*win, as.integer(end/win)*win, win)
  y <- unique(dlst[[4]])

  tab <- data.frame()
  
  nr <- length(x)
  nc <- length(y)
  
  if (length(idx.0) > 0) {
    m0 <- Matrix::sparseMatrix(
      i = match(dlst[[1]][idx.0], x),
      j = match(dlst[[4]][idx.0], y),
      x = dlst[[3]][idx.0],
      dims=c(nr,nc))

    rownames(m0) <- x
    colnames(m0) <- y
    m0 <- as.matrix(m0)
    tab1 <- reshape2::melt(m0)
    tab1$strand <- '+'
    tab <- rbind(tab,tab1)
  }

  if (length(idx.1) > 0) {
    m1 <- Matrix::sparseMatrix(
      i = match(dlst[[1]][idx.1], x),
      j = match(dlst[[4]][idx.1], y),
      x = dlst[[3]][idx.1],
      dims=c(nr,nc))

    rownames(m1) <- x
    colnames(m1) <- y
    m1 <- as.matrix(m1)
    tab2 <- reshape2::melt(m1)
    tab2$strand <- '-'
    tab <- rbind(tab,tab2)
  }

  colnames(tab) <- c("pos","label","depth","strand")
  tab$depth <- as.integer(tab$depth/win)
  return(tab)
}

#' @import dplyr
#' @export
plot.genes <- function(region = NULL, db = NULL, genes = NULL, label=TRUE,
                       max.overlaps = 5, collapse=TRUE)
{
  if (is.null(db)) stop(paste0("No database found, ", genome))
  
  genes1 <- db$gene
  exons1 <- db$exon
  trans1 <- db$transcript

  genes1 <- sort(sortSeqlevels(genes1), ignore.strand = TRUE)
  exons1 <- sort(sortSeqlevels(exons1), ignore.strand = TRUE)
  trans1 <- sort(sortSeqlevels(trans1), ignore.strand = TRUE)

  genes0 <- data.frame(subsetByOverlaps(genes1, region, ignore.strand = TRUE))
  genes <- genes %||% unique(genes0$gene_name)
  genes <- intersect(genes, genes0$gene_name)
  exons0 <- data.frame(subsetByOverlaps(exons1, region, ignore.strand = TRUE))
  exons0 <- exons0[which(exons0$gene_name %in% genes),]
  
  trans0 <- data.frame(subsetByOverlaps(trans1, region, ignore.strand = TRUE))
  trans0 <- trans0[which(trans0$gene_name %in% genes),]

  txcnt <- table(trans0$gene_name)
  trans0$idx = 1:nrow(trans0)
  tx = trans0$idx
  names(tx) = trans0$transcript_id

  exons0$idx = tx[exons0$transcript_id]

  trans0$start0 <- ifelse(trans0$strand == "-", trans0$end, trans0$start)
  trans0$end0 <- ifelse(trans0$strand == "-", trans0$start, trans0$end)

  require(ggarchery)
  p = ggplot()
  p = p + geom_arrowsegment(data = trans0,aes(x = start0, xend = end0, y = idx, yend = idx, color=strand),
                            arrow_positions = c(0.2, 0.4, 0.6, 0.8), 
                            arrows = arrow(length = unit(0.1, "inches")),
                            size = 1)
  p = p + geom_segment(data = exons0, aes(x = start, xend = end, y = idx, yend = idx, color=strand), size = 5)

  offset <- as.integer(max(trans0$idx)/5)
  offset <- ifelse(offset < 1, 1, offset)
  if (label & sum(as.character(trans0$strand) == '-') > 0) {
    p = p + ggrepel::geom_label_repel(
      data=trans0[which(as.character(trans0$strand)=="-"),], 
      aes(x = end, y = idx, label = gene_name), max.overlaps = max.overlaps,
      nudge_y = -trans0[which(as.character(trans0$strand)=="-"),]$idx-offset-1, size = 5, direction = "x")
  }
    
  
  if (label & sum(as.character(trans0$strand) == '+') > 0) {
    p = p + ggrepel::geom_label_repel(
      data=trans0[which(as.character(trans0$strand)=="+"),], 
      aes(x = start, y = idx, label = gene_name), max.overlaps = max.overlaps,
      nudge_y = -trans0[which(as.character(trans0$strand)=="+"),]$idx-offset, size = 5, direction = "x")
  }

  p = p + theme_void()  #+ facet_grid(gene_name~.) 
  p <- p + theme(panel.spacing= unit(0, "lines"),
                 axis.title.x=element_blank(), axis.text.y=element_blank(), 
                 axis.ticks.y=element_blank(), legend.position = "nona") +
    ylab("") + xlab("") +scale_x_continuous(limits=c(start(region), end(region)))

  offset <- -offset-2
  p + ylim(offset,max(trans0$idx+1)) + scale_color_manual(values = c("+" = "red", "-" = "blue"))
}

plot.bed <- function(region = NULL, peaks = NULL)
{
  r <- subsetByOverlaps(peaks, region, ignore.strand = TRUE)
  tab <- data.frame(r)
  
  p <- ggplot() + geom_segment(data = subset(tab, strand=="+"),aes(x = start, xend = end, y = 1, yend = 1, color=strand), size = 3)
  p <- p + geom_segment(data = subset(tab, strand=="-"),aes(x = start, xend = end, y = 0, yend = 0, color=strand), size = 3)
  p <- p + ylab("") + xlab("") + theme_void() + scale_x_continuous(limits=c(start(region), end(region)))
  p <- p + scale_color_manual(values = c("+" = "red", "-" = "blue"))
  return(p)
}
#' @import patchwork
#' @export
plot.cov <- function(bamfile=NULL, chr=NULL, start=-1, end =-1,
                     strand = c("both", "forward", "reverse", "ignore"),
                     max.depth = 0,
                     split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                     cell.group=NULL, log.scaled = log.scaled, start0 = -1, end0 = -1)
{
  if (is.null(chr) || start == -1 || end == -1) stop("Require a genomic region.")

  bc <- bamcov(bamfile=bamfile, chr=chr, start=start, end=end, strand=strand, split.bc=split.bc, cell.group=cell.group, bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)
  if (isTRUE(log.scaled)) {
    bc$depth <- log(bc$depth+1)
  }

  if (max.depth > 0) {
    bc$depth[bc$depth > max.depth] <- max.depth
  }
  
  bc$depth <- bc$depth * ifelse(bc$strand=='+',1,-1)

  ymax <- max(bc$depth)
  ymin <- min(bc$depth)
  
  p1 <- ggplot(bc, aes(x=pos,y=depth,fill=strand)) + geom_area(stat = "identity") +  facet_wrap(facets = ~label, strip.position = 'right', ncol = 1)  + ylab("") + theme_bw() + scale_x_continuous(limits=c(start, end))
  p1 <- p1 + scale_fill_manual(values = c("+" = "red", "-" = "blue")) + xlab("")
  
  if (start0 != -1 & start0 > start & end0 != -1 & end0 < end) {
    p1 <- p1 + annotate("rect", xmin = start0, xmax = end0, ymin = ymin, ymax = ymax, alpha = .1,fill = "blue")
  }
  return(p1)
}

#' @import patchwork
#' @export
plotTracks <-  function(bamfile=NULL, chr=NULL, start=-1, end =-1, gene=NULL,
                        strand = c("both", "forward", "reverse", "ignore"),
                        split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                        db = NULL, max.overlaps = 5, max.depth = 0,
                        cell.group=NULL, display.genes = NULL, toUCSC=FALSE, peaks =NULL,
                        log.scaled = FALSE, upstream = 1000, downstream = 1000, start0 = NULL,
                        end0 = NULL, anno.col = "blue", collapse = TRUE)
                        
{
  start0 <- start0 %||% start
  end0 <- end0 %||% end

  if (!is.null(gene)) {
    if (is.null(db)) stop(paste0("No database found, ", genome))
    genes1 <- db$gene
    start <- start(genes1[which(genes1$gene_name == gene)][1])
    end <- end(genes1[which(genes1$gene_name == gene)])

    if (start0 < start | end0 > end) stop("Conflict target region and gene.")
  }
  
  start <- start - upstream
  end <- end + downstream
  
  if (start0 < 0) start0 <- 1
  
  p1 <- plot.cov(bamfile=bamfile, chr=chr, start=start, end=end, strand=strand, split.bc=split.bc, bin=bin, cell.tag=cell.tag, umi.tag=umi.tag, cell.group=cell.group, log.scaled=log.scaled, start0 = start0, end0 = end0, max.depth = max.depth)
  p1 <- p1 + theme_pubr()# + theme(panel.spacing = unit(1, "lines"))
  
  if (!is.null(peaks)) {
    
    gr <- GRanges(seqnames=chr, ranges = IRanges(start = start, width = end-start))
    p2 <- plot.bed(region = gr, peaks = peaks)
    if (start0 > start & end0 < end) {
      p2 <- p2 + annotate("rect", xmin = start0, xmax = end0, ymin = 0, ymax = 1, alpha = .1,fill = anno.col)
    }
    
    if (toUCSC) chr <- paste0("chr",chr)
    gr <- GRanges(seqnames=chr, ranges = IRanges(start = start, width = end-start))
    p3 <- plot.genes(gr, db=db, genes=display.genes, max.overlaps = max.overlaps)
    p3 <- p3 + theme_pubr()
    if (start0 > start & end0 < end) {
      p3 <- p3 + annotate("rect", xmin = start0, xmax = end0, ymin = 0, ymax = 1, alpha = .1,fill = anno.col)
    }
    
    return(p2/p1/p3 + plot_layout(heights = c(1,8,2)))
  } 

  if (toUCSC) chr <- paste0("chr",chr)
  gr <- GRanges(seqnames=chr, ranges = IRanges(start = start, width = end-start))
  p2 <- plot.genes(gr, db=db, genes=display.genes, collapse=collapse)
  
  return(p1 / p2 + plot_layout(heights = c(8, 2)))  
}

#' @import RColorBrewer
#' @import pheatmap
#' @export
plotModHeatmap <- function(lc = NULL,
                           cluster_rows = FALSE,
                           cluster_cols =FALSE,
                           show_rownames = FALSE,
                           show_colnames = FALSE,
                           ...
                           )
{
  ann <- lc$module
  k <- length(unique(ann[['module']]))
  
  fig <- pheatmap(lc$LC,
                  annotation_row = ann,
                  cluster_rows = cluster_rows,
                  cluster_cols = cluster_cols,
                  show_rownames = show_rownames,
                  show_colnames = show_colnames,
                  ...)
  fig
}
