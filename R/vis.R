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

#'@importFrom gtools mixedsort
#'@importFrom ggrepel geom_label_repel
#'@importFrom viridis scale_color_viridis scale_fill_viridis
#'@import ggrepel
#'@import scales
#' 
FbtPlot0 <- function(tab = NULL,
                     col.by = NULL,
                     cols = NULL,
                     shape.by = NULL,
                     xlab = NULL, ylab = NULL,
                     point.label = NULL,  label.size=3,                     
                     zoom.in = FALSE,
                     gtf = NULL,
                     print.genes = NULL,
                     start = NULL,
                     end = NULL,
                     max.genes = 20,
                     layout_heights = c(3,2),
                     ...)
{
  if (isTRUE(zoom.in)) {
    if (nrow(tab) == 0) {
      p <- ggplot()
      if (is.null(gtf)) {
        return(p)
      } else {
        chr <- unique(tab$chr)
        if (is.null(start)) {
          start <- min(tab$start)
        }
        if (is.null(end) || end == -1) {
          end <- max(tab$start)
        }
        p1 <- plot.genes(chr = chr, start = start, end = end, gtf = gtf)
        return(p/p1 + plot_layout(heights=c(5,3)))
      }
    }

    chr <- unique(tab$chr)
    if (is.null(start)) {
      start <- min(tab$start)
    }
    if (is.null(end) || end == -1) {
      end <- max(tab$start)
    }
    message(paste0("Zoom in ", chr, " : ", start, "-", end))
    data <- tab
    data[["bp_cum"]] <- data[["start"]]
    xlab = NULL
  } else {
    data_cum <- tab %>% group_by(chr) %>% summarise(max_bp = max(start)) %>%
      mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% select(chr, bp_add)
    
    data <- tab %>% inner_join(data_cum, by = "chr") %>% mutate(bp_cum = start + bp_add)
    axis_set <- data %>% group_by(chr) %>% summarize(center = mean(bp_cum))
  }
  data$name <- tab$name
  
  data <- data %>% arrange(col.by)

  qval.min <- min(data$qval) - 1
  qval.max <- max(data$qval) + 1
  
  p <- ggplot(data) 
  if  (isFALSE(zoom.in)) {
    xi <- data_cum$bp_add
    xi <- xi[-1]
    p <- p + geom_vline(xintercept = xi, color="red", linetype="dotted")
  }
  if (!is.null(col.by)) {
    if (is.null(shape.by)) {
      p <- p + geom_point(aes(x=bp_cum, y=qval, fill=.data[[col.by]]), shape=21, ...)
      if (is.numeric(data[[col.by]])) {
        p <- p + scale_fill_viridis()
      } else {
        p <- p + scale_fill_manual(values = cols)
      }
    } else {
      p <- p + geom_point(aes(x=bp_cum, y=qval, fill=.data[[col.by]], shape=.data[[shape.by]]), ...)
      p <- p + scale_shape_manual(values = c(21, 24, 22, 23, 25)) 
      if (is.numeric(data[[col.by]])) {
        p <- p + scale_fill_viridis()
      } else {
        p <- p + scale_fill_manual(values = cols)
      }
    }
  } else {
    if (is.null(shape.by)) {
      p <- p + geom_point(aes(x=bp_cum, y=qval),  ...)
    } else {
      p <- p + geom_point(aes(x=bp_cum, y=qval, shape = .data[[shape.by]]),  ...)
    }
  }
  if (isFALSE(zoom.in)) {
    p <- p + scale_x_continuous(label = axis_set$chr, breaks = axis_set$center,
                                limits = c(min(data$bp_cum), max(data$bp_cum)),
                                expand=c(0.01,0.01)) # , guide = guide_axis(n.dodge=2))
  } else {
    p <- p + coord_cartesian(xlim=c(start, end), ylim=c(qval.min, qval.max), expand=FALSE) + scale_x_continuous(labels = scales::label_comma())
  }
  p <- p + fbt_theme() + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))
  p <- p + xlab(xlab) + ylab(ylab)

  if (!is.null(point.label)) {
    sel <- intersect(point.label, data$name)
    if (length(sel) > 0) {
      p <- p + geom_label_repel(data=subset(data, name %in% sel),aes(x=bp_cum, y=qval,label=name),box.padding = 0.5, max.overlaps = Inf, size=label.size)
    }
  }
  if (isTRUE(zoom.in) && !is.null(gtf)) {
    chr <- unique(tab$chr)
    p1 <- plot.genes(chr = chr, start = start, end = end, gtf = gtf, print.genes = print.genes, max.genes = max.genes)
    return(p/p1 + plot_layout(heights=layout_heights))
  }
  p
}

#'@export
FbtPlot <- function(object = NULL,
                    val = NULL,
                    assay = NULL,
                    chr.name = "chr", start.name = "start", end.name = "end",
                    col.by = NULL, cols = NULL, sel.chrs = NULL,
                    xlab = "Chromosome", ylab = expression(-log[10](p)),
                    types = NULL,
                    point.label = NULL,
                    label.size=3,
                    idents = NULL,
                    shape.by= NULL,
                    chr = NULL, start = NULL, end = NULL,
                    gtf = NULL, gene = NULL, upstream=1000, downstream=1000,
                    print.genes = NULL,
                    layout_heights = c(3,2),
                    ...)
{
  if (is.null(val)) stop("No value name specified.")  
  cols <- cols %||% c("#131313","blue","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
  assay <- assay %||% DefaultAssay(object)

  n <- length(assay)

  chr1 <- chr
  start1 <- start
  end1 <- end
  
  zoom.in <- FALSE
  if (!is.null(chr)) {
    zoom.in <- TRUE
    if (is.null(start)) {
      start1 <- 0
    }
    if (is.null(end)) {
      end1 <- -1
    }
  }

  if (!is.null(gtf)) {
    if (!isGTF(gtf)) stop("Not like a GTF database, use gtf2db to read GTF first.")
    if (!is.null(gene)) {
      sl <- .Call("G_query_gene", gtf, gene)
      if (is.null(sl)) stop(paste0("No gene ", gene, " found in the database!!"))
      chr1 <- sl[[1]]
      start1 <- sl[[2]]
      end1 <- sl[[3]]
      zoom.in <- TRUE
    }
  }
  
  sl <- lapply(1:n, function(i) {
    assay0 <- assay[i]
    object0 <- object[[assay0]]
    tab0 <- object0[[]]
    
    if (chr.name %ni% colnames(tab0)) stop("No chr name found.")
    if (start.name %ni% colnames(tab0)) stop("No start name found.")
    if (val %ni% colnames(tab0)) stop("No val name found.")

    if (zoom.in) {
      tab0 %>% filter(chr %in% chr1 & start >= start1) -> tab0
      if (end1 > 0) {
        tab0 <- subset(tab0, start <=end1)
      }
    }

    if (!is.null(sel.chrs)) {
      tab0 <- tab0 %>% filter (chr %in% sel.chrs)
      levels(tab0[[chr.name]]) <- levels(sel.chrs)
    }
    
    if (!is.null(types)) {
      if ("type" %ni% colnames(tab0)) stop("No type found.")
      tab0 <- subset(tab0, type %in% types)
      if (nrow(tab0) == 0) stop("Empty records.")
    }
    if (is.null(levels(tab0[[chr.name]]))) {
      lv <- mixedsort(unique(tab0[[chr.name]]))
    } else {
      lv <- levels(tab0[[chr.name]])
    }
    tab <- data.frame(chr = factor(tab0[[chr.name]], levels = lv),
                      start = as.numeric(tab0[[start.name]]),
                      qval = -log10(as.numeric(tab0[[val]])),
                      name = rownames(tab0))

    if (end.name %in% colnames(tab0)) {
      tab[['start']] <- as.numeric(as.integer((tab0[[start.name]] + tab0[[end.name]])/2))
    }
    if (!is.null(col.by)) {
      tab[[col.by]] <- tab0[[col.by]]
    }

    if (!is.null(shape.by)) {
      tab[[shape.by]] <- tab0[[shape.by]]
    }

    tab[['assay']] <- assay0

    tab <- subset(tab, !is.na(qval))

    tab
  })

  tab <- data.table::rbindlist(sl)

  if (isTRUE(zoom.in)) {
    if (start1 > 0) {
      start1 <- start1 - upstream
    }
    if (end1 != -1) {
      end1 <- end1 + downstream
    }
  }
  if (n == 1) {
    p <- FbtPlot0(tab=tab, col.by=col.by, cols=cols, xlab=xlab, ylab = ylab, point.label=point.label, shape.by=shape.by, label.size=label.size, zoom.in = zoom.in, start = start1, end = end1, gtf = gtf, print.genes = print.genes, layout_heights = layout_heights, ...)
  } else {
    if (is.null(shape.by)) {
      shape.by <- "assay"
    }
    p <- FbtPlot0(tab=tab, col.by=col.by, cols = cols, shape.by = shape.by, xlab=xlab, ylab = ylab, point.label=point.label, label.size=label.size, zoom.in = zoom.in, start = start1, end = end1, gtf = gtf, print.genes = print.genes, layout_heights = layout_heights, ...)
  }
  p
}

theme_cov <- function(...) {
  theme(
    axis.title.y = element_text(color = "black", family = "Helvetica",size = rel(1)),
    axis.title.x = element_blank(),
    axis.text.y = element_text(family = "Helvetica",color = "black",size = rel(1)),
    axis.text.x = element_text(family = "Helvetica",color = "black",size = rel(1.5)),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "whitesmoke"),
    plot.title = element_blank(),
    ...
  )
}

plot.genes <- function(chr = NULL, start = NULL, end = NULL, gtf = NULL, genes = NULL, label=TRUE, highlights=NULL, print.genes = NULL, max.genes = 20, ...)
{
  if (is.null(gtf)) stop("No database specified.")
  if (!isGTF(gtf)) stop("Not like a GTF database, use gtf2db to read GTF first.")
  if (is.null(chr)) stop("No chromosome specified.")
  
  if (is.null(start)) {
    start = 0;
  }
  if (is.null(end)) {
    end = -1;
  }
  
  tracks <- .Call("gene_tracks", as.character(chr), start, end, gtf, genes)
  
  if (is.null(tracks)) {
    p <- ggplot()
    p <- p + ylab("") + xlab("")
    if (end > 0) {
      p <- p + coord_cartesian(xlim=c(start, end), expand=FALSE)
    }
    p <- p + ylim(0,1)
    
    if (!is.null(highlights)) {
      df <- as.data.frame(highlights)
      df$ymin <- 0
      df$ymax <- 1
      p <- p + geom_rect(data=df, inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax,ymin=ymin,ymax=ymax), color="grey", alpha=0.2)
    }
    p <- p + fbt_theme()
    return(p)
  }

  idx <- which(tracks$start < start)
  tracks[idx,"start"] <- start
  idx <- which(tracks$end >end)
  if (end > 0) {
    tracks[idx,"end"] <- end
    tracks <- subset(tracks, start < end)
  }

  mi <- as.integer(max(tracks$idx)/3 *4)
  tracks %>%
    filter(type == 1) %>%
    group_by(gene) %>%
    summarise(start = min(start), end = max(end), idx = max(idx), nudge_y = mi - max(idx)) -> gname
  
  gname <- as.data.frame(gname)
  gname$med <- (gname$start + gname$end)/2

  if (!is.null(print.genes)) {
    gname %>% filter(gene %in% print.genes)-> gname
  } else {
    if (nrow(gname) > max.genes) {
      if (max.genes > 0) {
        message(paste0("Over ", max.genes, " genes. Only print ", max.genes, " gene names."))
        idx <- sample(1:nrow(gname), max.genes)
        gname <- gname[idx,]
      }
    }
  }
  p <- ggplot() + geom_segment(data = subset(tracks, type == 1),aes(x = start, xend = end, y = idx, yend = idx, color=strand), size=1)
  p <- p + geom_segment(data = subset(tracks, type == 2), aes(x = start, xend = end, y = idx, yend = idx),color="black", size = 3)

  if (!is.null(highlights)) {
    df <- as.data.frame(highlights)
    df$ymin <- 0
    df$ymax <- mi
    p <- p + geom_rect(data=df,inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax,ymin=ymin,ymax=ymax), color="grey", alpha=0.2)
  }
  if (max.genes > 0 || !is.null(print.genes)) {
    p <- p + geom_text_repel(data=gname,aes(x=med,y=idx, label=gene), nudge_y = gname$nudge_y, size=5, max.overlaps=Inf, segment.color = "grey50")
  }
  p <- p + fbt_theme()
  p <- p + theme(panel.spacing= unit(0, "lines"), axis.text = element_blank(),
                 axis.title =element_blank(), 
                 axis.ticks =element_blank())
  p <- p + ylab("") + xlab("") + coord_cartesian(xlim=c(start, end), expand=FALSE)
  p <- p + ylim(0,mi)
  p <- p + scale_color_manual(values = c("+" = "red", "-" = "blue"))
  p
}

plot.bed <- function(start = NULL, end = NULL, peaks = NULL, type.col = NULL, group.title.size=rel(2), highlights=NULL)
{
  tab <- subset(peaks, start >= start, end <= end)
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
  p <- p + xlab("") + coord_cartesian(xlim=c(start, end),expand=FALSE)
  p <- p + ylab("") + scale_y_continuous(expand = c(0,0))
  p <- p + theme_void() + theme(legend.position = "none",
                                panel.background = element_rect(fill = "grey"),
                                panel.spacing= unit(0.1, "lines"),
                                strip.text = element_text(size = group.title.size))
  p <- p + ylab("EPTs")

  return(p)
}

#' @import patchwork
#' @importFrom scales pretty_breaks
plot.cov <- function(bamfile=NULL, chr=NULL, start=-1, end =-1,
                     strand = c("both", "forward", "reverse", "ignore"),
                     max.depth = 0,
                     split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                     cell.group=NULL, log.scaled = log.scaled, #start0 = -1, end0 = -1,
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
plot.cov2 <- function(fragfile=NULL, chr=NULL, start=-1, end =-1,
                      max.depth = 0,
                      split.bc = FALSE, bin = 1000,
                      cell.group=NULL, log.scaled = FALSE) #, start0 = -1, end0 = -1)
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
#' @export
plotTracks <-  function(bamfile=NULL, chr=NULL, start=NULL, end =NULL, gene=NULL,
                        strand = "both",
                        split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                        gtf = NULL, max.depth = 0, group.title.size = rel(2),
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
    if (is.null(gtf)) stop("gtf is not specified.")
    if (!isGTF(gtf)) stop("Not like a GTF database, use gtf2db to read GTF first.")
    sl <- .Call("G_query_gene", gtf, gene)
    if (is.null(sl)) stop(paste0("No gene ", gene, " found in the database!!"))
    
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
    p0 <- plot.bed(start = start, end = end, peaks = tab, type.col=type.col, highlights = df)
  } 

  p3 <- NULL

  if (!is.null(fragfile)) {
    p3 <- plot.cov2(fragfile = fragfile, chr=chr, start=start, end=end, bin=bin, cell.group=cell.group, log.scaled = atac.log.scaled, max.depth=atac.max.depth)
    p3 <- p3 + theme_cov()
    p3 <- p3 + theme(panel.spacing.y = unit(0.1, "lines"))
    p3 <- p3 + theme(strip.text = element_text(size = group.title.size))
  }
  
  p2 <- plot.genes(chr = chr, start = start, end = end, gtf=gtf, genes=display.genes, collapse=collapse, highlights=df, ...)
  
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
