#'@importFrom dplyr %>%
#'@importFrom gtools mixedsort
#'@importFrom ggrepel geom_label_repel
#'@importFrom viridis scale_color_viridis scale_fill_viridis
#'@import ggrepel
#' 
#'@export
FbtPlot0 <- function(tab = NULL, col.by = NULL, cols = NULL, shape.by = NULL, xlab = "Chromosome", ylab = expression(-log[10](P)), point.label = NULL, arrange.type = FALSE, label.size=3, ...)
{
  data_cum <- tab %>% group_by(chr) %>% summarise(max_bp = max(start)) %>%
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% select(chr, bp_add)
  
  data <- tab %>% inner_join(data_cum, by = "chr") %>% mutate(bp_cum = start + bp_add)
  axis_set <- data %>% group_by(chr) %>% summarize(center = mean(bp_cum))

  data$name <- tab$name
  #rownames(data) <- rownames(tab)
  
  if (isTRUE(arrange.type)) data <- data %>% arrange(type)

  fbt_theme <- function() {
    theme(
      legend.text = element_text(face = "italic",color = "black",family = "Helvetica",size = rel(1.5)),
      axis.title.y = element_text(color = "black", family = "Helvetica",size = rel(1)),
      axis.title.x = element_text(color = "black", family = "Helvetica",size = rel(1.5)),
      axis.text = element_text(family = "Helvetica",color = "black",size = rel(1.5)),
      axis.line = element_blank(),
      axis.ticks = element_line(size = rel(1), color = "black"),
      panel.border = element_rect(color = "black", fill = NA, size= rel(2), linetype = "solid"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "whitesmoke"),
      legend.key = element_rect(fill = "whitesmoke"),
      legend.title = element_text(size = rel(1.5),family = "Helvetica"),
      plot.title = element_text(color = "black",face = "bold",size = rel(1.7),family = "Helvetica")
    )
  }

  xi <- data_cum$bp_add
  xi <- xi[-1]
  p <- ggplot(data) + geom_vline(xintercept = xi, color="red", linetype="dotted")
  
  if (!is.null(col.by)) {
    if (is.null(shape.by)) {
      p <- p + geom_point(aes(x=bp_cum, y=pval, fill=.data[[col.by]]), shape=21, ...)
      if (is.numeric(data[[col.by]])) {
        p <- p + scale_fill_viridis()
      } else {
        p <- p + scale_fill_manual(values = cols)
      }
    } else {
      p <- p + geom_point(aes(x=bp_cum, y=pval, col=.data[[col.by]], shape=.data[[shape.by]]), ...)
      if (is.numeric(data[[col.by]])) {
        p <- p + scale_color_viridis()
      } else {
        p <- p + scale_color_manual(values = cols)
      }
    }
  } else {
    if (is.null(shape.by)) {
      p <- p + geom_point(aes(x=bp_cum, y=pval),  ...)
    } else {
      p <- p + geom_point(aes(x=bp_cum, y=pval, shape = .data[[shape.by]]),  ...)
    }
  }
  p <- p + scale_x_continuous(label = axis_set$chr, breaks = axis_set$center,
                              limits = c(min(data$bp_cum), max(data$bp_cum)),
                              expand=c(0.01,0.01))
  p <- p + fbt_theme() + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))
  p <- p + xlab(xlab) + ylab(ylab)

  if (!is.null(point.label)) {
    sel <- intersect(point.label, data$name)
    if (length(sel) > 0) {
      p <- p + geom_label_repel(data=subset(data, name %in% sel),aes(x=bp_cum, y=pval,label=name),box.padding = 0.5, max.overlaps = Inf, size=label.size)
    }
  }
  p
}

#'@export
FbtPlot <- function(object = NULL, assay = NULL, chr = "chr", start = "start", val = NULL, col.by = NULL, cols = NULL, sel.chrs = NULL, xlab = "Chromosome", ylab = expression(-log[10](P)), types = NULL, point.label = NULL, arrange.type = FALSE, label.size=3, idents = NULL, ...)
{
  if (is.null(val)) stop("No value name specified.")  
  cols <- cols %||% c("#131313","blue","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
  assay <- assay %||% DefaultAssay(object)

  n <- length(assay)

  sl <- lapply(1:n, function(i) {
    assay0 <- assay[i]
    object0 <- object[[assay0]]
    tab0 <- object0[[]]
    
    if (chr %ni% colnames(tab0)) stop("No chr name found.")
    if (start %ni% colnames(tab0)) stop("No start name found.")
    if (val %ni% colnames(tab0)) stop("No val name found.")
    
    if (!is.null(sel.chrs)) {
      tab0 <- tab0 %>% filter (chr %in% sel.chrs)
    }
    
    if (!is.null(types)) {
      if ("type" %ni% colnames(tab0)) stop("No type found.")
      tab0 <- subset(tab0, type %in% types)
      if (nrow(tab0) == 0) stop("Empty records.")
    }
    lv <- mixedsort(unique(tab0[[chr]]))
    
    tab <- data.frame(chr = factor(tab0[[chr]], levels = lv),
                      start = as.numeric(tab0[[start]]),
                      pval = -log10(as.numeric(tab0[[val]])),
                      name = rownames(tab0))
    
    if (!is.null(col.by)) {
      tab[[col.by]] <- tab0[[col.by]]
    }
    
    tab[['assay']] <- assay0

    tab <- subset(tab, !is.na(pval))
    tab
  })

  tab <- data.table::rbindlist(sl)

  if (n == 1) {
    p <- FbtPlot0(tab=tab, col.by=col.by, cols=cols, xlab=xlab, ylab = ylab, point.label=point.label, arrange.type = FALSE, label.size=label.size, ...)
  } else {
    p <- FbtPlot0(tab=tab, col.by=col.by, cols = cols, shape.by = "assay", xlab=xlab, ylab = ylab, point.label=point.label, arrange.type = FALSE, label.size=label.size,  ...)
  }
  p
}

theme_cov <- function(...) {
  theme(
    #legend.text = element_blank(),
    axis.title.y = element_text(color = "black", family = "Helvetica",size = rel(1)),
    axis.title.x = element_blank(),
    axis.text.y = element_text(family = "Helvetica",color = "black",size = rel(1)),
    axis.text.x = element_text(family = "Helvetica",color = "black",size = rel(1.5)),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "whitesmoke"),
    #legend.position = "none",
    plot.title = element_blank(),
    ...
  )
}

#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomeInfoDb sortSeqlevels
#' @export
plot.genes <- function(region = NULL, db = NULL, genes = NULL, label=TRUE, highlights=NULL, ...)
{
  if (is.null(db)) stop("No database specified.")
  if (!isGTF(db)) stop("Not like a GTF database, use gtf2db to read GTF first.")
  
  start0 <- start(region)
  end0 <- end(region)

  tracks <- .Call("gene_tracks", as.character(seqnames(region)), start0, end0, db, genes)
  
  if (is.null(tracks)) {
    p <- ggplot()
    p <- p + ylab("") + xlab("") + coord_cartesian(xlim=c(start(region), end(region)), expand=FALSE)
    p <- p + ylim(0,1) + theme_void()
    
    if (!is.null(highlights)) {
      df <- as.data.frame(highlights)
      df$ymin <- 0
      df$ymax <- 1
      p <- p + geom_rect(data=df, inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax,ymin=ymin,ymax=ymax), color="grey", alpha=0.2)
    }
    return(p)
  }

  idx <- which(tracks$start < start0)
  tracks[idx,"start"] <- start0
  idx <- which(tracks$end >end0)
  tracks[idx,"end"] <- end0
  tracks <- subset(tracks, start < end)

  mi <- as.integer(max(tracks$idx)/3 *4)
  tracks %>%
    filter(type == 1) %>%
    group_by(gene) %>%
    summarise(start = min(start), end = max(end), idx = max(idx), nudge_y = mi - max(idx)) -> gname
  
  gname <- as.data.frame(gname)
  gname$med <- (gname$start + gname$end)/2
  
  p <- ggplot() + geom_segment(data = subset(tracks, type == 1),aes(x = start, xend = end, y = idx, yend = idx, color=strand), size=1)
  p <- p + geom_segment(data = subset(tracks, type == 2), aes(x = start, xend = end, y = idx, yend = idx),color="black", size = 3)

  if (!is.null(highlights)) {
    df <- as.data.frame(highlights)
    df$ymin <- 0
    df$ymax <- mi
    p <- p + geom_rect(data=df,inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax,ymin=ymin,ymax=ymax), color="grey", alpha=0.2)
  }

  p <- p + geom_text_repel(data=gname,aes(x=med,y=idx, label=gene), nudge_y = gname$nudge_y, size=5, max.overlaps=Inf, segment.color = "grey50")
  p <- p + theme_minimal()
  p <- p + theme(panel.spacing= unit(0, "lines"), axis.text = element_blank(),
                 axis.title =element_blank(), 
                 axis.ticks =element_blank())
  p <- p + ylab("") + xlab("") + coord_cartesian(xlim=c(start(region), end(region)), expand=FALSE)
  p <- p + ylim(0,mi)
  p <- p + scale_color_manual(values = c("+" = "red", "-" = "blue"))
  p
}
#'@importFrom IRanges subsetByOverlaps
plot.bed <- function(region = NULL, peaks = NULL, type.col = NULL, group.title.size=rel(2), highlights=NULL)
{
  r <- subsetByOverlaps(peaks, region, ignore.strand = TRUE)
  tab <- data.frame(r)
  p <- ggplot()
  if ("type" %in% colnames(tab) & !is.null(type.col)) {
    p <- p + geom_rect(data = tab,inherit.aes = F,aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill=type), size = 1)
    #p <- p + geom_segment(data = subset(tab, strand=="-"),aes(x = start, xend = end, y = 0, yend = 0, color=type), size = 3)
    p <- p + scale_fill_manual(values = type.col)
  } else {
    p <- p + geom_rect(data = tab, inherit.aes = F,aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill=strand), size = 1)
    #p <- p + geom_segment(data = subset(tab, strand=="-"),aes(x = start, xend = end, y = 0, yend = 0, color=strand), size = 3)
    p <- p + scale_fill_manual(values = c("+" = "red", "-" = "blue"))
  }
  if (!is.null(highlights)) {
    df <- as.data.frame(highlights)
    df$ymin <- 0
    df$ymax <- 1
    p <- p + geom_rect(data=df,inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax,ymin=ymin,ymax=ymax), color="grey", alpha=0.2)
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
#' @importFrom GenomicRanges seqnames
#' @importFrom scales pretty_breaks
#' @export
plot.cov <- function(bamfile=NULL, chr=NULL, start=-1, end =-1,
                     strand = c("both", "forward", "reverse", "ignore"),
                     max.depth = 0,
                     split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                     cell.group=NULL, log.scaled = log.scaled, start0 = -1, end0 = -1,
                     highlights=NULL,
                     junc=FALSE, junc.min.depth = 10)
{
  if (is.null(chr) || start == -1 || end == -1) stop("Require a genomic region.")
  
  bc <- bamcov(bamfile=bamfile, chr=as.character(chr), start=start, end=end, strand=strand,
               split.bc=split.bc,
               cell.group=cell.group, bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)

  if (junc) {
    juncs <- bamjunc(bamfile=bamfile, chr=as.character(chr), start=start, end=end, strand=strand,
                     split.bc=split.bc,
                     cell.group=cell.group, cell.tag=cell.tag, umi.tag=umi.tag)
  }
  
  if (!is.null(cell.group)) {
    ss <- table(unlist(cell.group))
    bc$depth <- bc$depth/as.vector(ss[as.character(bc$label)])
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

  posmax <- max(bc$pos)
  posmin <- min(bc$pos)

  p1 <- ggplot() + geom_area(data=bc, aes(x=pos,y=depth,fill=strand), stat = "identity")
  
  if (junc) {
    juncs <- subset(juncs, depth >= junc.min.depth)
        
    if (!is.null(cell.group)) {
      juncs$depth <- juncs$depth / as.vector(ss[as.character(juncs$label)])
    }

    if (isTRUE(log.scaled)) {
      juncs$depth <- log1p(juncs$depth)
    }

    ymax0 <- max(abs(bc$depth))
    juncs$depth <- juncs$depth/ymax0
    juncs$depth[which(juncs$depth>1)] <- 1
    juncs[["y"]] <- 0

    juncs <- subset(juncs, end > start)
    idx <- which (juncs$strand == "-")
    
    juncs["depth"][idx,] <- juncs["depth"][idx,] * -1
    
    p1 <- p1 + geom_splice(data=juncs, aes(x=start, xend = end, y = y, height = depth))
  }
  
  if (!is.null(highlights)) {
    df <- as.data.frame(highlights)
    df$ymin <- ymin
    df$ymax <- ymax
    p1 <- p1 + geom_rect(data=df,inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax,ymin=ymin,ymax=ymax), color="grey", alpha=0.2)
  }
  
  p1 <- p1 + scale_y_continuous(breaks=pretty_breaks(),guide=guide_axis(check.overlap = T))
  p1 <- p1 + facet_wrap(facets = ~label, strip.position = 'right', ncol = 1)
  p1 <- p1 + xlab("") + ylab("") + theme_bw() +coord_cartesian(xlim=c(start, end), expand=FALSE)
  p1 <- p1 + scale_fill_manual(values = c("+" = "red", "-" = "blue"))
  p1 <- p1 + theme_cov()
  p1 <- p1 + theme(panel.spacing.y = unit(0.1, "lines"))
  p1 <- p1 + theme(legend.position="none")

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
  bc <- fragcov(fragfile=fragfile, chr=as.character(chr), start=start, end=end,
                split.bc=split.bc, cell.group=cell.group, bin=bin)

  if (!is.null(cell.group)) {
    ss <- table(unlist(cell.group))
    bc$depth <- bc$depth/as.vector(ss[as.character(bc$label)])
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
                        anno.col = "blue", type.col = NULL, layout_heights =c(1,10,2),
                        highlights = NULL, junc = FALSE,
                        ...)
                        
{
  if (!is.null(gene)) {
    if (is.null(db)) stop("db is not specified.")
    if (!isGTF(db)) stop("Not like a GTF database, use gtf2db to read GTF first.")
    sl <- .Call("G_query_gene", db, gene)

    chr <- sl[[1]]
    start <- sl[[2]]
    end <- sl[[3]]
  }

  if (is.null(start) || is.null(end)) stop("No start or/and end position specified.")
  start <- start - upstream
  end <- end + downstream
  message(paste0("Chromosome ", chr, ", start ", start, ", end ", end))
  if (is.null(bamfile)) stop("No bam file specified.")

  df <- NULL
  if (!is.null(highlights)) {
    if (is.list(highlights)) {
      highlights <- lapply(highlights, function(x) x[c(1,2)])
    }

    df <- matrix(unlist(highlights),nrow=2)
    df <- as.data.frame(t(df))
    colnames(df) <- c("xmin","xmax")

    df <- subset(df, xmin > start & xmax < end)
  }
  
  p1 <- plot.cov(bamfile=bamfile, chr=chr, start=start, end=end, strand=strand, split.bc=split.bc,
                 bin=bin, cell.tag=cell.tag, umi.tag=umi.tag, cell.group=cell.group,
                 log.scaled=log.scaled, max.depth = max.depth, highlights=df, junc=junc)
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
    p0 <- plot.bed(region = gr, peaks = peaks, type.col=type.col, highlights = df)

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
  p2 <- plot.genes(gr, db=db, genes=display.genes, collapse=collapse, highlights=df, ...)
  
  if (!is.null(p0)) {
    if (!is.null(p3)) {
      return(p0/ p1 / p3/ p2 + plot_layout(heights=layout_heights[c(1,2,2,3)]))
    } else {
      return(p0/ p1 / p2 + plot_layout(heights=layout_heights))
    }
  }
  if (!is.null(p3)) {
    return(p1 / p3/ p2 + plot_layout(heights=layout_heights[c(2,2,3)]))
  }
  return(p1 / p2 + plot_layout(heights=layout_heights[c(2,3)]))
}

## #' @import RColorBrewer
## #' @importFrom pheatmap pheatmap
## #' @export
## plotModHeatmap <- function(lc = NULL,
##                            cluster_rows = FALSE,
##                            cluster_cols =FALSE,
##                            show_rownames = FALSE,
##                            show_colnames = FALSE,
##                            ...
##                            )
## {
##   ann <- lc$module
##   k <- length(unique(ann[['module']]))  
##   fig <- pheatmap(lc$LC,
##                   annotation_row = ann,
##                   cluster_rows = cluster_rows,
##                   cluster_cols = cluster_cols,
##                   show_rownames = show_rownames,
##                   show_colnames = show_colnames,
##                   ...)
##   fig
## }
