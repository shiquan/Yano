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
  message(paste0("Process ", bamfile))

  strand.flag <- -1
  if (strand == "forward") strand.flag <- 0
  else if(strand == "reverse") strand.flag <- 1
  else if (strand == "ignore") strand.flag <- -2

  if (!is.null(cell.group)) {
    #if (length(cell.group)) stop("cell.group should be a vector of cell groups with cell names.")    
    cell.names <- names(cell.group)
    if (is.null(cell.names)) stop("Names of cell.group should not be empty.")
    cell.group <- as.character(cell.group)
    
    groups <- levels(cell.group) %||% unique(cell.group)
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
#'@useDynLib Yano depth2matrix
#'@export
fragcov <- function(fragfile = NULL, chr = NULL, start = -1, end = -1, split.bc = FALSE, cell.group = NULL, bin=1000)
{
  if (is.null(fragfile)) stop("No fragment file.")
  if (is.null(chr)) stop("No chromosome name")
  if (start == -1) stop("No start position.")
  if (end == -1) stop("No end position")

  len <- end - start
  if (len > 10000000) { ## 10M
    stop("Too large region, try to shrink it before querying.")
  }

  win <- 1
  if (len > bin*2) {
    win <- as.integer(len/bin)
  }
  message(paste0("Process ", fragfile))

  if (!is.null(cell.group)) {
    #if (length(cell.group)) stop("cell.group should be a vector of cell groups with cell names.")    
    cell.names <- names(cell.group)
    if (is.null(cell.names)) stop("Names of cell.group should not be empty.")
    cell.group <- as.character(cell.group)
    
    groups <- levels(cell.group) %||% unique(cell.group)
    group.ids <- match(cell.group, groups)
    n <- length(cell.group)
    split.bc <- TRUE    
    dlst <- .Call("fragment2matrix", fragfile, chr, start, end, cell.names, n, group.ids, groups)
  } else {
    dlst <- .Call("fragment2matrix", fragfile, chr, start, end,  NULL, 0, NULL, NULL)
  }
  
  dlst[[1]] <- as.integer(dlst[[1]]/win)*win

  x <- seq(as.integer(start/win)*win, as.integer(end/win)*win, win)
  y <- unique(dlst[[4]])

  nr <- length(x)
  nc <- length(y)
  
  m0 <- Matrix::sparseMatrix(
    i = match(dlst[[1]], x),
    j = match(dlst[[4]], y),
    x = dlst[[3]],
    dims=c(nr,nc))

  rownames(m0) <- x
  colnames(m0) <- y
  m0 <- as.matrix(m0)
  tab <- reshape2::melt(m0)
  tab$strand <- '.'

  colnames(tab) <- c("pos","label","depth","strand")
  tab$depth <- as.integer(tab$depth/win)

  return(tab)
}

theme_cov <- function(...) {
  theme(
    legend.text = element_blank(),
    axis.title.y = element_text(color = "black", family = "Helvetica",size = rel(1)),
    axis.title.x = element_blank(),
    axis.text = element_text(family = "Helvetica",color = "black",size = rel(1.5)),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "whitesmoke"),
    legend.position = "none",
    plot.title = element_blank(),
    ...
  )
}

#' @import dplyr
#' @importFrom IRanges subsetByOverlaps
#' @export
plot.genes <- function(region = NULL, db = NULL, genes = NULL, label=TRUE, ...)
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
  if (nrow(trans0) == 0) {
    p <- ggplot()
    p <- p + ylab("") + xlab("") + coord_cartesian(xlim=c(start(region), end(region)), expand=FALSE)
    #scale_x_continuous(limits=c(start(region), end(region)),expand=c(0,0))
    p <- p + ylim(0,1) + theme_void()
    return(p)
  }

  trans0$idx <- .Call("trans_sort",trans0$start, trans0$end)
  tx = trans0$idx
  names(tx) = trans0$transcript_id
  exons0$idx = tx[exons0$transcript_id]

  start0 <- start(region)
  end0 <- end(region)
  trans0$start[which(trans0$start<start0)] <- start0
  trans0$end[which(trans0$end>end0)] <- end0
  gname <- trans0 %>% group_by(gene_name) %>% summarise(start = min(start), end = max(end), idx = max(idx)+1)
  gname <- as.data.frame(gname)
  gname$med <- (gname$start + gname$end)/2
  p <- ggplot() + geom_segment(data = trans0,aes(x = start, xend = end, y = idx, yend = idx, color=strand), size=1)
  p <- p + geom_segment(data = exons0, aes(x = start, xend = end, y = idx, yend = idx),color="black", size = 5)
  p <- p + geom_text(data=gname,aes(x=med,y=idx,label=gene_name), size=5, check_overlap = TRUE,na.rm=TRUE)
  p <- p + theme_minimal()
  p <- p + theme(panel.spacing= unit(0, "lines"), axis.text = element_blank(),
                 axis.title =element_blank(), 
                 axis.ticks =element_blank(),
                 legend.position="left",
                 legend.justification="right",
                 legend.box.spacing = unit(-10, "pt"),
                 legend.margin=margin(0,0,0,0))
  p <- p + ylab("") + xlab("") + coord_cartesian(xlim=c(start(region), end(region)), expand=FALSE)
  # scale_x_continuous(limits=c(start(region), end(region)), expand = c(0,0))
  p <- p + ylim(0,max(trans0$idx)+1)
  p <- p + scale_color_manual(values = c("+" = "red", "-" = "blue"))
  p 
}
#'@importFrom IRanges subsetByOverlaps
plot.bed <- function(region = NULL, peaks = NULL, type.col = NULL, group.title.size=rel(2))
{
  r <- subsetByOverlaps(peaks, region, ignore.strand = TRUE)
  tab <- data.frame(r)
  p <- ggplot()
  if ("type" %in% colnames(tab) & !is.null(type.col)) {
    p <- p + geom_segment(data = tab,aes(x = start, xend = end, y = 0, yend = 1, color=type), size = 1)
    #p <- p + geom_segment(data = subset(tab, strand=="-"),aes(x = start, xend = end, y = 0, yend = 0, color=type), size = 3)
    p <- p + scale_color_manual(values = type.col)
  } else {
    p <- p + geom_segment(data = tab, aes(x = start, xend = end, y = 0, yend = 1, color=strand), size = 1)
    #p <- p + geom_segment(data = subset(tab, strand=="-"),aes(x = start, xend = end, y = 0, yend = 0, color=strand), size = 3)
    p <- p + scale_color_manual(values = c("+" = "red", "-" = "blue"))
  }
  p <- p + facet_wrap(facets = ~strand, strip.position = 'right', ncol = 1)
  p <- p + xlab("") + coord_cartesian(xlim=c(start(region), end(region)),expand=FALSE)
  #scale_x_continuous(limits=c(start(region), end(region)), expand = c(0,0))
  p <- p + ylab("") + scale_y_continuous(expand = c(0,0))
  p <- p + theme_void() + theme(legend.position = "none",
                                panel.background = element_rect(fill = "grey"),
                                panel.spacing= unit(0.1, "lines"),
                                strip.text = element_text(size = group.title.size))
  p <- p + ylab("EPTs")

  return(p)
}
#' @import patchwork
#' @import dplyr
#' @importFrom GenomicRanges seqnames
#' @export
plot.cov <- function(bamfile=NULL, chr=NULL, start=-1, end =-1,
                     strand = c("both", "forward", "reverse", "ignore"),
                     max.depth = 0,
                     split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                     cell.group=NULL, log.scaled = log.scaled, start0 = -1, end0 = -1)
{
  if (is.null(chr) || start == -1 || end == -1) stop("Require a genomic region.")

  if (is.list(bamfile)) {
    nm <- names(bamfile)
    if (is.list(cell.group)) {
      dl <- lapply(nm, function(x) {
        bamcov(bamfile=bamfile[[x]], chr=as.character(chr), start=start, end=end, strand=strand,
               split.bc=split.bc,
               cell.group=cell.group[[x]], bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)
      })
    } else {
      dl <- lapply(nm, function(x) {
        bamcov(bamfile=bamfile[[x]], chr=as.character(chr), start=start, end=end, strand=strand,
               split.bc=split.bc,
               cell.group=cell.group, bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)
      })
    }

    bc <- bind_rows(dl) %>%  group_by(pos, label, strand) %>% summarise(sum(depth, na.rm = TRUE))
    ss <- table(unlist(bc))
    colnames(bc) <- c("pos", "label", "strand", "depth")
  } else {
    bc <- bamcov(bamfile=bamfile, chr=as.character(chr), start=start, end=end, strand=strand, split.bc=split.bc,
                 cell.group=cell.group, bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)
  }

  bc$label <- as.character(bc$label)
  bc$label <- factor(bc$label, levels=gtools::mixedsort(unique(bc$label)))
  if (!is.null(cell.group)) {
    ss <- table(unlist(cell.group))
    bc$depth <- bc$depth/as.vector(ss[bc$label]) *1000
  }
  
  if (isTRUE(log.scaled)) {
    bc$depth <- log1p(bc$depth)
  }

  if (max.depth > 0) {
    bc$depth[bc$depth > max.depth] <- max.depth
  }
  
  bc$depth <- bc$depth * ifelse(bc$strand=='+',1,-1)

  ymax <- max(bc$depth)
  ymin <- min(bc$depth)
  
  p1 <- ggplot(bc, aes(x=pos,y=depth,fill=strand)) + geom_area(stat = "identity")
  p1 <- p1 + facet_wrap(facets = ~label, strip.position = 'right', ncol = 1)
  p1 <- p1 + xlab("") + ylab("") + theme_bw() +coord_cartesian(xlim=c(start, end), expand=FALSE)
  # scale_x_continuous(limits=c(start, end),expand=c(0,0))
  p1 <- p1 + scale_fill_manual(values = c("+" = "red", "-" = "blue"))  
  ## if (start0 != -1 & start0 > start & end0 != -1 & end0 < end) {
  ##   p1 <- p1 + annotate("rect",
  ##                       xmin = start0,
  ##                       xmax = end0,
  ##                       ymin = ymin,
  ##                       ymax = ymax,
  ##                       alpha = .1,
  ##                       fill = "blue")
  ## }
  return(p1)
}

#' @import patchwork
#' @import dplyr
#' @importFrom GenomicRanges seqnames
#' @export
plot.cov2 <- function(fragfile=NULL, chr=NULL, start=-1, end =-1,
                      max.depth = 0,
                      split.bc = FALSE, bin = 1000,
                      cell.group=NULL, log.scaled = FALSE, start0 = -1, end0 = -1)
{
  if (is.null(chr) || start == -1 || end == -1) stop("Require a genomic region.")

  if (is.list(fragfile)) {
    nm <- names(fragfile)
    if (is.list(cell.group)) {
      dl <- lapply(nm, function(x) {
        fragcov(fragfile=fragfile[[x]], chr=as.character(chr), start=start, end=end,
                split.bc=split.bc,cell.group=cell.group[[x]], bin=bin)
      })
    } else {
      dl <- lapply(nm, function(x) {
        fragcov(fragfile=fragfile[[x]], chr=as.character(chr), start=start, end=end,
                split.bc=split.bc, cell.group=cell.group, bin=bin)
      })
    }

    bc <- bind_rows(dl) %>%  group_by(pos, label, strand) %>% summarise(sum(depth, na.rm = TRUE))
    ss <- table(unlist(bc))
    colnames(bc) <- c("pos", "label", "strand", "depth")
  } else {
    bc <- fragcov(fragfile=fragfile, chr=as.character(chr), start=start, end=end,
                  split.bc=split.bc, cell.group=cell.group, bin=bin)
  }
  
  bc$label <- as.character(bc$label)
  bc$label <- factor(bc$label, levels=gtools::mixedsort(unique(bc$label)))
  if (!is.null(cell.group)) {
    ss <- table(unlist(cell.group))
    bc$depth <- bc$depth/as.vector(ss[bc$label]) *1000
  }
  
  if (isTRUE(log.scaled)) {
    bc$depth <- log1p(bc$depth)
  }

  if (max.depth > 0) {
    bc$depth[bc$depth > max.depth] <- max.depth
  }
  
  ymax <- max(bc$depth)
  ymin <- min(bc$depth)
  
  p1 <- ggplot(bc, aes(x=pos,y=depth), fill="black") + geom_area(stat = "identity")
  p1 <- p1 + facet_wrap(facets = ~label, strip.position = 'right', ncol = 1)
  p1 <- p1 + xlab("") + ylab("") + theme_bw() +coord_cartesian(xlim=c(start, end), expand=FALSE)

  return(p1)
}

#' @import patchwork
#' @importFrom GenomicRanges seqnames GRanges start end
#' @export
plotTracks <-  function(bamfile=NULL, chr=NULL, start=NULL, end =NULL, gene=NULL,
                        strand = "both",
                        split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                        db = NULL, max.depth = 0, group.title.size = rel(2),
                        cell.group=NULL, display.genes = NULL, toUCSC=FALSE, meta.features =NULL,
                        log.scaled = FALSE, upstream = 1000, downstream = 1000,
                        fragfile = NULL,
                        atac.log.scaled = FALSE,
                        atac.max.depth = 0,
                        anno.col = "blue", type.col = NULL, layout_heights =c(1,10,10,2),...)
                        
{
  if (!is.null(gene)) {
    if (is.null(db)) stop(paste0("No database found, ", genome))
    genes1 <- db$gene

    start <- start %||% min(start(genes1[which(genes1$gene_name %in% gene)]))
    end <- end %||% max(end(genes1[which(genes1$gene_name %in% gene)]))
    chr <- unique(seqnames(genes1[which(genes1$gene_name %in% gene)]))
    if (length(chr) != 1) stop(paste("More than 1 chromosome found, ", chr))
  }
  if (is.null(start) || is.null(end)) stop("No start or/and end position specified.")
  start <- start - upstream
  end <- end + downstream
  message(paste0("chr ", chr, ", start ", start, ", end ", end))
  if (is.null(bamfile)) stop("No bam file specified.")
    
  p1 <- plot.cov(bamfile=bamfile, chr=chr, start=start, end=end, strand=strand, split.bc=split.bc,
                 bin=bin, cell.tag=cell.tag, umi.tag=umi.tag, cell.group=cell.group,
                 log.scaled=log.scaled, max.depth = max.depth)
  p1 <- p1 + theme_cov()
  p1 <- p1 + theme(panel.spacing.y = unit(0.1, "lines"))
  p1 <- p1 + theme(strip.text = element_text(size = group.title.size))

  p0 <- NULL
  if (!is.null(meta.features)) {
    if (length(intersect(c("chr","start","end","strand","type"), colnames(meta.features))) != 5)
      stop("No chr/start/end/strand/type column found in meta.features")

    tab <- subset(meta.features, start > 0 & end > 0)
    peaks <- GRanges(seqnames=tab[['chr']],
                 ranges = IRanges(start = tab[['start']],
                                  end = tab[['end']]),
                 strand = tab[['strand']],
                 type = tab[['type']])

    gr <- GRanges(seqnames=chr, ranges = IRanges(start = start, width = end-start))
    p0 <- plot.bed(region = gr, peaks = peaks, type.col=type.col)

    ## if (start0 > start & end0 < end) {
    ##   p2 <- p2 + annotate("rect", xmin = start0, xmax = end0, ymin = 0, ymax = 1, alpha = .1,fill = anno.col)
    ## }    
  } 

  p3 <- NULL

  if (!is.null(fragfile)) {
    p3 <- plot.cov2(fragfile = fragfile, chr=chr, start=start, end=end, bin=bin, cell.group=cell.group, log.scaled = atac.log.scaled, max.depth=atac.max.depth)
    p3 <- p3 + theme_cov()
    p3 <- p3 + theme(panel.spacing.y = unit(0.1, "lines"))
    p3 <- p3 + theme(strip.text = element_text(size = group.title.size))
  }
  
  if (toUCSC) chr <- paste0("chr",chr)
  gr <- GRanges(seqnames=chr, ranges = IRanges(start = start, width = end-start))
  p2 <- plot.genes(gr, db=db, genes=display.genes, collapse=collapse, ...)
  
  if (!is.null(p0)) {
    if (!is.null(p3)) {
      return(p0/ p1 / p3/ p2 + plot_layout(heights=layout_heights))
    } else {
      return(p0/ p1 / p2 + plot_layout(heights=layout_heights[c(1,2,4)]))
    }
  }
  if (!is.null(p3)) {
    return(p1 / p3/ p2 + plot_layout(heights=layout_heights[c(2,3,4)]))
  }
  return(p1 / p2 + plot_layout(heights=layout_heights[c(2,4)]))
}

#' @import RColorBrewer
#' @importFrom pheatmap pheatmap
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
