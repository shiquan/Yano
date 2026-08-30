# Safe font family: fall back to system default if Helvetica is not available
.fbt_family <- function() {
  if (nzchar(system.file(package = "systemfonts"))) {
    avail <- systemfonts::system_fonts()$family
    if ("Helvetica" %in% avail) return("Helvetica")
  }
  ""  # ggplot2 default sans-serif
}

fbt_theme <- function() {
  family <- .fbt_family()
  theme(
    legend.text = element_text(face = "italic",color = "black",family = family,size = rel(1.5)),
    axis.title.y = element_text(color = "black", family = family,size = rel(1)),
    axis.title.x = element_text(color = "black", family = family,size = rel(1.5)),
    axis.text = element_text(family = family,color = "black",size = rel(1.5)),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = rel(1), color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = rel(2), linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "whitesmoke"),
    legend.key = element_rect(fill = "whitesmoke"),
    legend.title = element_text(size = rel(1.5),family = family),
    plot.title = element_text(color = "black",face = "bold",size = rel(1.7),family = family)
  )
}

#'@importFrom gtools mixedsort
#'@importFrom ggrepel geom_label_repel
#'@importFrom viridis scale_color_viridis scale_fill_viridis
#'@import ggrepel
#'@importFrom scales label_comma

# Internal helper: add discrete or continuous fill scale
.add_fill_scale <- function(p, data, col.by, cols) {
  if (is.numeric(data[[col.by]])) {
    return(p + scale_fill_viridis())
  }
  n <- length(unique(data[[col.by]]))
  if (n < 15) {
    cols <- cols %||% c("#131313","blue","#A6CEE3","#1F78B4","#B2DF8A",
                        "#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00",
                        "#CAB2D6","#6A3D9A","#FFFF99","#B15928")
  } else {
    cols <- cols %||% colours(distinct = TRUE)[seq_len(n)]
  }
  p + scale_fill_manual(values = cols)
}

FbtPlot0 <- function(tab = NULL,
                     col.by = NULL,
                     cols = NULL,
                     shape.by = NULL,
                     xlab = NULL, ylab = NULL,
                     point.label = NULL,  label.size=3,
                     pt.size = NULL,
                     zoom.in = FALSE,
                     gtf = NULL,
                     print.genes = NULL,
                     start = NULL,
                     end = NULL,
                     max.genes = 20,
                     layout.heights = c(3,2),                     
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
        return(p/p1 + plot_layout(heights=layout.heights))
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
    data_cum <- tab %>% group_by(chr) %>% summarise(max_bp = max(start), .groups = "drop") %>%
      mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% select(chr, bp_add)
    
    data <- tab %>% inner_join(data_cum, by = "chr") %>% mutate(bp_cum = start + bp_add)
    axis_set <- data %>% group_by(chr) %>% summarize(center = mean(bp_cum))
  }
  data$name <- tab$name
  
  if (!is.null(col.by)) {
    data <- data %>% arrange(.data[[col.by]])
  }

  qval.min <- min(data$qval) - 1
  qval.max <- max(data$qval) + 1
  if (qval.min < -0.01 & qval.max < 10) {
    qval.min <- -0.01
  }


  pt.size <- pt.size %||% 1
  p <- ggplot(data) 
  if  (isFALSE(zoom.in)) {
    xi <- data_cum$bp_add
    xi <- xi[-1]
    p <- p + geom_vline(xintercept = xi, color="grey60", linetype="dotted")
  }
  if (!is.null(col.by)) {
    if (is.null(shape.by)) {
      p <- p + geom_point(aes(x=bp_cum, y=qval, fill=.data[[col.by]]), shape=21, size = pt.size)
      p <- .add_fill_scale(p, data, col.by, cols)
    } else {
      p <- p + geom_point(aes(x=bp_cum, y=qval, fill=.data[[col.by]], shape=.data[[shape.by]]), size = pt.size)
      p <- p + scale_shape_manual(values = c(21, 24, 22, 23, 25))
      p <- .add_fill_scale(p, data, col.by, cols)
      p <- p + guides(fill = guide_legend(override.aes = list(shape=21)))
    }
  } else {
    if (is.null(shape.by)) {
      p <- p + geom_point(aes(x=bp_cum, y=qval),  size = pt.size)
    } else {
      p <- p + geom_point(aes(x=bp_cum, y=qval, shape = .data[[shape.by]]),  size = pt.size)      
    }
  }
  x_axis_angle <- 0
  if (isFALSE(zoom.in)) {
    ## Assess chromosome label overlap risk using width distribution
    chr_widths <- tab %>%
      group_by(chr) %>%
      summarise(width = max(start), .groups = "drop")
    total_width <- sum(chr_widths$width)
    if (total_width == 0) total_width <- 1  # degenerate guard
    chr_widths$width_frac <- chr_widths$width / total_width
    chr_widths$name_len <- nchar(as.character(chr_widths$chr))
    ## squeeze: chars per unit of relative width (higher = more overlap risk)
    chr_widths$squeeze <- chr_widths$name_len / chr_widths$width_frac
    max_squeeze <- max(chr_widths$squeeze)
    mean_squeeze <- sum(chr_widths$name_len)
    squeeze_ratio <- max_squeeze / mean_squeeze
    if (is.infinite(squeeze_ratio) || is.nan(squeeze_ratio)) squeeze_ratio <- 999

    n_chrs <- nrow(axis_set)
    max_chr_len <- max(chr_widths$name_len)

    if (squeeze_ratio < 3 && n_chrs <= 12 && max_chr_len <= 4) {
      x_axis_text_size <- rel(1.5)
      x_axis_guide <- waiver()
    } else if (squeeze_ratio < 10 && n_chrs <= 22) {
      x_axis_text_size <- rel(1.2)
      x_axis_guide <- waiver()
    } else if (squeeze_ratio < 40) {
      x_axis_text_size <- rel(1.0)
      x_axis_guide <- guide_axis(n.dodge = 2)
    } else {
      x_axis_text_size <- rel(0.85)
      x_axis_guide <- waiver()
      x_axis_angle <- 45
    }
    p <- p + scale_x_continuous(labels = axis_set$chr, breaks = axis_set$center,
                                limits = c(min(data$bp_cum), max(data$bp_cum)),
                                expand = c(0.01, 0.01),
                                guide = x_axis_guide)
  } else {
    x_axis_text_size <- rel(1.5)
    p <- p + coord_cartesian(xlim=c(start, end), ylim=c(qval.min, qval.max), expand=FALSE) + scale_x_continuous(labels = scales::label_comma())
  }
  p <- p + fbt_theme() + theme(axis.title.y = element_text(size = rel(1.5), angle = 90))
  p <- p + theme(axis.text.x = element_text(size = x_axis_text_size,
                                             angle = x_axis_angle,
                                             hjust = if (x_axis_angle > 0) 1 else 0.5))
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
    p <- p + theme(panel.spacing= unit(0, "lines"), axis.text.x = element_blank(), axis.title.x =element_blank())
    return(p/p1 + plot_layout(heights=layout.heights))
  }
  p
}
#' @title FbtPlot
#' @description Generate Manhatten plot for log10 scaled p value for spatial dissimilarity test. X axis is the location of features. The function use ggplot2 as backend.
#' @param object Seurat object.
#' @param val Specify the name of p value. The name usually be format like bind-name.pval and bind-name.padj (BH adjusted p value). For example, if you use gene as binding feature, the name should be gene_name.pval or gene_name.padj.
#' @param assay Work assay.
#' @param chr.name The title of chromosome name in the meta table. Default is "chr".
#' @param start.name The title of start position name in the meta table. Default is "start".
#' @param end.name The title of end position name in the meta table. Default is "end".
#' @param color.by Color points by specify the title of values in meta table. Can be discrete or continous.
#' @param col.by Alias to "color.by".
#' @param colors Manually specify the colors. Used with color.by.
#' @param cols Alias to "colors".
#' @param sel.chrs Vector of selected chromosome names to plot. Change the order by set the level of chr names.
#' @param pt.size Point size.
#' @param xlab Label for x axis. Default is "Chromosome".
#' @param ylab Label for y axis. Default is "-log10p".
#' @param subset Rule for subsetting the meta table before plot.
#' @param point.label Vector of points to plot their labels.
#' @param label.size Size of label. Default is 3.
#' @param shape.by Shape points by specify the title of values in meta table. Can only be discrete.
#' @param chr Choromsome to zoom in. Default is NULL, no zoom in. The zoom in mode can be enabled by setting chr or gene/gtf.
#' @param start Start position to zoom in.
#' @param end End position to zoom in.
#' @param gtf GTF database. Load by \code{gtf2db}. Default is NULL. If specified transcirpt tracks will be plotted.
#' @param gene Gene name to zoom in. Should used with gtf database specified.
#' @param upstream Flank zoom in region with upstream. Default is 1000. Only works when zoom in mode enabled.
#' @param downstream Flank zoom in region with downstream. Default is 1000. Only works when zoom in mode enabled.
#' @param print.genes Print the gene names in the transcript tracks. Default will print all or randomly 20 genes if more than 20 genes in this region.
#' @param remove.chr Remove 'chr' in the chromosome names.
#' @param layout.heights Specify the layouts for Manhatten plot and gene tracks. Default is c(3,2).
# #' @param ... Parameters pass to geom_point().
#' 
#'@export
FbtPlot <- function(object = NULL,
                    val = NULL,
                    assay = NULL,
                    chr.name = "chr", start.name = "start", end.name = "end",
                    color.by = NULL, colors = NULL,
                    col.by = NULL, cols = NULL,
                    sel.chrs = NULL,
                    pt.size = NULL,
                    xlab = "Chromosome", ylab = expression(-log[10](p)),
                    subset = NULL,
                    point.label = NULL,
                    label.size=3,
                    shape.by= NULL,
                    chr = NULL, start = NULL, end = NULL,
                    gtf = NULL, gene = NULL, upstream=1000, downstream=1000,
                    print.genes = NULL,
                    remove.chr = FALSE,
                    layout.heights = c(3,2))
{
  if (is.null(val)) stop("No value name specified.")

  col.by <- col.by %||% color.by
  cols <- cols %||% colors
  
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

    if (isTRUE(zoom.in)) {
      tab0 <- base::subset(tab0, tab0[[chr.name]] %in% chr1 & tab0[[start.name]] >= start1)
      if (end1 > 0) {
        tab0 <- base::subset(tab0, tab0[[start.name]] <= end1)
      }
    }

    if (isTRUE(remove.chr)) {
      tab0[[chr.name]] <- gsub("^chr", "", as.character(tab0[[chr.name]]))
      if (!is.null(sel.chrs)) {
        sel.chrs <- gsub("^chr", "", as.character(sel.chrs))
      }
    }
    if (!is.null(sel.chrs)) {
      tab0[[chr.name]] <- as.character(tab0[[chr.name]])
      sel.chrs <- as.character(sel.chrs)
      tab0 <- tab0 %>% filter(.data[[chr.name]] %in% sel.chrs)
      tab0[[chr.name]] <- factor(tab0[[chr.name]], levels = sel.chrs)
    }
    
    ## if (!is.null(types)) {
    ##   if ("type" %ni% colnames(tab0)) stop("No type found.")
    ##   tab0 <- subset(tab0, type %in% types)
    ##   if (nrow(tab0) == 0) stop("Empty records.")
    ## }
    if (!is.null(subset)) {
      tab0 <- base::subset(tab0, subset = subset)
      if (nrow(tab0) == 0) stop("Empty records after subset.")
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
      tab[['start']] <- as.numeric((tab0[[start.name]] + tab0[[end.name]]) / 2)
    }
    tab[['assay']] <- assay0

    if (!is.null(col.by)) {
      if (col.by %in% colnames(tab0)) {
        tab[[col.by]] <- tab0[[col.by]]
      } else if (!col.by %in% colnames(tab)) {
        warning(paste0("No ", col.by, " found at meta.features."))
      }
    }

    if (!is.null(shape.by)) {
      if (shape.by %in% colnames(tab0)) {
        tab[[shape.by]] <- tab0[[shape.by]]
      } else if (!shape.by %in% colnames(tab)) {
        warning(paste0("No ", shape.by, " found at meta.features."))
      }
    }
    
    tab <- base::subset(tab, !is.na(qval))

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
  p <- NULL
  if (n == 1) {
    p <- FbtPlot0(tab=tab, col.by=col.by, cols=cols, xlab=xlab, ylab = ylab, point.label=point.label, shape.by=shape.by, label.size=label.size, zoom.in = zoom.in, start = start1, end = end1, gtf = gtf, print.genes = print.genes, layout.heights = layout.heights, pt.size = pt.size)
  } else {
    p <- FbtPlot0(tab=tab, col.by=col.by, cols = cols, shape.by = shape.by, xlab=xlab, ylab = ylab, point.label=point.label, label.size=label.size, zoom.in = zoom.in, start = start1, end = end1, gtf = gtf, print.genes = print.genes, layout.heights = layout.heights, pt.size = pt.size)
  }
  p
}

theme_cov <- function(...) {
  family <- .fbt_family()
  theme(
    axis.title.y = element_text(color = "black", family = family, size = rel(1)),
    axis.title.x = element_blank(),
    axis.text.y = element_text(family = family, color = "black", size = rel(0.8)),
    axis.text.x = element_text(family = family, color = "black", size = rel(1.5)),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_rect(fill = "whitesmoke"),
    plot.title = element_blank(),
    ...
  )
}
#' @export
plot.genes <- function(x = NULL, chr = NULL, start = NULL, end = NULL, gtf = NULL, genes = NULL, label=TRUE, highlights=NULL, print.genes = NULL, max.genes = 20, ...)
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
      p <- p + geom_rect(data=df, inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax), ymin=-Inf, ymax=Inf, fill="grey", alpha=0.2)
      p <- p + geom_vline(data=df, inherit.aes = F, mapping=aes(xintercept=xmin), color="darkorange", linewidth=0.6, alpha=0.9)
    }
    p <- p + fbt_theme()
    return(p)
  }

  idx <- which(tracks$start < start)
  tracks[idx, "start"] <- start
  if (end > 0) {
    idx <- which(tracks$end > end)
    tracks[idx, "end"] <- end
  }
  tracks <- subset(tracks, start < end)

  if (nrow(tracks) == 0) {
    return(ggplot() + theme_void())
  }

  mi <- as.integer(max(tracks$idx)/3 *4)
  if (mi < max(tracks$idx)+1) {
    mi <-  max(tracks$idx)+1
  }
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
  p <- ggplot() + geom_segment(data = subset(tracks, type == 1),aes(x = start, xend = end, y = idx, yend = idx, color=strand), linewidth=1)
  p <- p + geom_segment(data = subset(tracks, type == 2), aes(x = start, xend = end, y = idx, yend = idx),color="black", linewidth = 3)

  if (!is.null(highlights)) {
    df <- as.data.frame(highlights)
    p <- p + geom_rect(data=df, inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax), ymin=-Inf, ymax=Inf, fill="grey", alpha=0.2)
    p <- p + geom_vline(data=df, inherit.aes = F, mapping=aes(xintercept=xmin), color="darkorange", linewidth=0.6, alpha=0.9)
  }
  if (max.genes > 0 || !is.null(print.genes)) {
    p <- p + geom_text_repel(data=gname,aes(x=med,y=idx, label=gene), nudge_y = gname$nudge_y, size=5, max.overlaps=Inf, segment.color = "grey50")
  }
  p <- p + scale_x_continuous(labels = scales::label_comma())
  p <- p + theme(panel.spacing= unit(0, "lines"), axis.text.y = element_blank(), axis.title.y =element_blank(), axis.ticks.y =element_blank())
  p <- p + ylab("") + xlab("") + coord_cartesian(xlim=c(start, end), expand=FALSE)
  p <- p + ylim(0,mi)
  p <- p + scale_color_manual(values = c("+" = "red", "-" = "blue", "." = "grey60"))
  p <- p + fbt_theme()
  p
}
#' @export
plot.bed <- function(x = NULL, start = NULL, end = NULL, peaks = NULL, col.by = NULL, group.title.size=rel(2), highlights=NULL, ...)
{
  start.filter <- start
  end.filter <- end
  tab <- subset(peaks, start >= start.filter & end <= end.filter)
  p <- ggplot()
  if ("type" %in% colnames(tab) & !is.null(col.by)) {
    p <- p + geom_rect(data = tab,inherit.aes = F,aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill=type), linewidth = 1)
    #p <- p + geom_segment(data = subset(tab, strand=="-"),aes(x = start, xend = end, y = 0, yend = 0, color=type), linewidth = 3)
    p <- p + scale_fill_manual(values = col.by)
  } else {
    p <- p + geom_rect(data = tab, inherit.aes = F,aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill=strand), linewidth = 1)
    #p <- p + geom_segment(data = subset(tab, strand=="-"),aes(x = start, xend = end, y = 0, yend = 0, color=strand), linewidth = 3)
    p <- p + scale_fill_manual(values = c("+" = "red", "-" = "blue", "." = "grey60"))
  }
  if (!is.null(highlights)) {
    df <- as.data.frame(highlights)
    p <- p + geom_rect(data=df, inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax), ymin=-Inf, ymax=Inf, fill="grey", alpha=0.2)
    p <- p + geom_vline(data=df, inherit.aes = F, mapping=aes(xintercept=xmin), color="darkorange", linewidth=0.6, alpha=0.9)
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

#' @importFrom scales pretty_breaks
plot.cov <- function(bamfile=NULL, chr=NULL, start=-1, end =-1,
                     strand = c("both", "forward", "reverse", "ignore"),
                     max.depth = 0,
                     split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                     cell.group=NULL, log.scaled = FALSE, #start0 = -1, end0 = -1,
                     highlights=NULL,
                     junc=FALSE, junc.min.depth = 0)
{
  if (is.null(chr) || start == -1 || end == -1) stop("Require a genomic region.")
  if (is.list(cell.group)) {
    nm <- names(cell.group)
    if (!is.null(nm)) {
      if (is.list(bamfile)) {
        nm2 <- names(bamfile)
        if (!is.null(nm2)) {
          bad.nm <- setdiff(nm, nm2)
          if (length(bad.nm) > 0) {
            warning("Inconsistance list name between cell group and bam files.")
          }
          nm <- intersect(nm, nm2)
          bamfile <- bamfile[nm]
        } 
      }
    }
  }
  bc <- bamcov(bamfile=bamfile, chr=as.character(chr), start=start, end=end, strand=strand,
               split.bc=split.bc,
               cell.group=cell.group, bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)
  if (is.null(bc)) {
    message("No reads found in region ", chr, ":", start, "-", end, ".")
    return(ggplot())
  }
  if (junc) {
    juncs <- bamjunc(bamfile=bamfile, chr=as.character(chr), start=start, end=end, strand=strand,
                     split.bc=split.bc,
                     cell.group=cell.group, cell.tag=cell.tag, umi.tag=umi.tag)
    if (is.null(juncs)) {
      juncs <- data.frame(start = integer(0), end = integer(0),
                          label = character(0), depth = numeric(0),
                          strand = character(0))
    }
  }

  if (!is.null(cell.group)) {
    ss <- table(unlist(cell.group))
    ss0 <- ss[as.character(bc$label)]
    ss0[is.na(ss0)] <- 1  # group absent from table (unused factor level)
    bc$depth <- bc$depth/as.vector(ss0)
  }
  
  if (isTRUE(log.scaled)) {
    bc$depth <- log1p(bc$depth)
  }

  bc$depth <- bc$depth * ifelse(bc$strand=='-',-1,1)
  
  ymax <- max(bc$depth)
  ymin <- min(bc$depth)
  ymin0 <- ymin
  if (ymin < 0) {
    ymin0 <- ifelse(max.depth > 0, -1*max.depth, ymin)
  }
  ymax0 <- ifelse(max.depth > 0, max.depth, ymax)
  y_pad <- (ymax0 - ymin0) * 0.06
  ymin0 <- ymin0 - y_pad
  ymax0 <- ymax0 + y_pad
  
  posmax <- max(bc$pos)
  posmin <- min(bc$pos)

  if (strand == "ignore") {
    bc$depth <- abs(bc$depth)
    # Collapse the two per-strand records into one summed value per
    # (pos, label). Keeping both would make geom_area draw a
    # self-intersecting polygon and leave alternating unfilled bands.
    bc <- bc %>% group_by(pos, label) %>% summarise(depth = sum(depth), .groups = "drop")
    bc["strand"] <- "."
  }

  p1 <- ggplot() + geom_area(data=bc, aes(x=pos,y=depth,fill=strand), stat = "identity")
  
  if (junc) {

    juncs <- subset(juncs, depth >= junc.min.depth)
    if (nrow(juncs) > 0) {
      if (!is.null(cell.group)) {
        juncs$depth <- juncs$depth / as.vector(ss[as.character(juncs$label)])
      }
      
      if (isTRUE(log.scaled)) {
        juncs$depth <- log1p(juncs$depth)
      }

      juncs[["y"]] <- 0
      
      juncs <- subset(juncs, end > start)
      juncs <- subset(juncs, depth > 0)

      idx <- which (juncs$strand == "-")
      
      juncs$depth[idx] <- juncs$depth[idx] * -1
      
      p1 <- p1 + geom_splice(data=juncs, aes(x=start, xend = end, y = y, height = depth))
    }
  }
  
  if (!is.null(highlights)) {
    df <- as.data.frame(highlights)
    p1 <- p1 + geom_rect(data=df, inherit.aes = F, mapping=aes(xmin=xmin, xmax=xmax), ymin=-Inf, ymax=Inf, fill="grey", alpha=0.2)
    p1 <- p1 + geom_vline(data=df, inherit.aes = F, mapping=aes(xintercept=xmin), color="darkorange", linewidth=0.6, alpha=0.9)
  }

  p1 <- p1 + scale_y_continuous(breaks = pretty_breaks(n = 4), labels = label_comma(), guide = guide_axis(check.overlap = TRUE), expand = expansion(mult = c(0.06, 0.06)))
  p1 <- p1 + facet_wrap(facets = ~label, strip.position = 'right', ncol = 1)
  p1 <- p1 + xlab("") + ylab("") + theme_bw() +coord_cartesian(xlim=c(start, end), ylim = c(ymin0, ymax0), expand=FALSE)
  p1 <- p1 + scale_fill_manual(values = c("+" = "red", "-" = "blue", "." = "grey60"))
  p1 <- p1 + theme_cov()
  p1 <- p1 + theme(panel.spacing.y = unit(0.1, "lines"))
  p1 <- p1 + theme(legend.position="none")

  return(p1)
}

#' @import dplyr
plot.cov2 <- function(fragfile=NULL, chr=NULL, start=-1, end =-1,
                      max.depth = 0,
                      split.bc = FALSE, bin = 1000,
                      cell.group=NULL, log.scaled = FALSE) #, start0 = -1, end0 = -1)
{
  if (is.null(chr) || start == -1 || end == -1) stop("Require a genomic region.")
  bc <- fragcov(fragfile=fragfile, chr=as.character(chr), start=start, end=end,
                split.bc=split.bc, cell.group=cell.group, bin=bin)

  if (is.null(bc)) {
    message("No fragments found in region ", chr, ":", start, "-", end, ".")
    return(ggplot())
  }

  if (!is.null(cell.group)) {
    ss <- table(unlist(cell.group))
    ss0 <- ss[as.character(bc$label)]
    ss0[is.na(ss0)] <- 1  # group absent from table (unused factor level)
    bc$depth <- bc$depth/as.vector(ss0)
  }
  
  if (isTRUE(log.scaled)) {
    bc$depth <- log1p(bc$depth)
  }

  if (max.depth > 0) {
    bc$depth[bc$depth > max.depth] <- max.depth
  }
  
  ymax <- max(bc$depth)
  ymin <- min(bc$depth)
  
  p1 <- ggplot(bc, aes(x=pos,y=depth)) + geom_area(stat = "identity", fill="black")
  p1 <- p1 + facet_wrap(facets = ~label, strip.position = 'right', ncol = 1)
  p1 <- p1 + scale_y_continuous(expand = expansion(mult = c(0.06, 0.06)))
  p1 <- p1 + xlab("") + ylab("") + theme_bw() +coord_cartesian(xlim=c(start, end), expand=FALSE)

  return(p1)
}

# Add cross-panel highlight rectangles and boundary lines to a gtable using grid primitives.
#
# After patchwork assembles sub-plots into a combined gtable, this function finds all
# panel grobs and adds spanning rectGrob + segmentsGrob for each highlight region.
# This bypasses the per-panel squish_infinite() transform that limits
# geom_rect(ymin=-Inf, ymax=Inf) to a single facet panel.
#
# Handles two panel types produced by patchwork:
#   1. gtable  – faceted sub-plot (e.g., plot.cov ~label): the inner gtable has
#                sub-panels (panel-1-1, panel-1-2, ...).  Span them all.
#   2. gTree   – non-faceted sub-plot (e.g., plot.genes): add highlight grobs as
#                children of the panel gTree.
#' @noRd
.add_grid_highlights_to_gtable <- function(gt, highlights, start, end) {
  layout <- gt$layout

  # Find ALL panel entries at the top level (across all sub-plots)
  panel_idx <- which(grepl("^panel", layout$name) & !grepl("panel-area", layout$name))
  if (length(panel_idx) == 0L || nrow(highlights) == 0L) return(gt)

  # Combined bounding box of ALL panels from ALL sub-plots
  t_all <- min(layout$t[panel_idx])
  b_all <- max(layout$b[panel_idx])
  l_all <- min(layout$l[panel_idx])
  r_all <- max(layout$r[panel_idx])

  # Protrude beyond panels for visual prominence
  extend_rows <- 2L
  t_span <- max(1L, t_all - extend_rows)
  b_span <- b_all + extend_rows

  window <- end - start

  # Build highlight rect + boundary-line grobs once
  y_extend <- 0.03  # 3% overhang on each side
  hl_grobs <- vector("list", nrow(highlights) * 2L)
  for (i in seq_len(nrow(highlights))) {
    x0 <- (highlights$xmin[i] - start) / window
    x1 <- (highlights$xmax[i] - start) / window

    hl_grobs[[2L * i - 1L]] <- rectGrob(
      x      = unit(x0, "npc"),
      width  = unit(x1 - x0, "npc"),
      y      = unit(-y_extend, "npc"),
      height = unit(1 + 2 * y_extend, "npc"),
      just   = c("left", "bottom"),
      gp     = gpar(fill = "grey", alpha = 0.2, col = NA)
    )
    hl_grobs[[2L * i]] <- segmentsGrob(
      x0 = unit(x0, "npc"), x1 = unit(x0, "npc"),
      y0 = unit(-y_extend, "npc"), y1 = unit(1 + y_extend, "npc"),
      gp = gpar(col = "darkorange", lwd = 0.6, alpha = 0.9)
    )
  }
  hl_tree <- do.call(grobTree, hl_grobs)

  # Single spanning grob across ALL panels from ALL tracks
  gt <- gtable::gtable_add_grob(gt, hl_tree,
    t = t_span, l = l_all, b = b_span, r = r_all,
    clip = "off", name = "hl_span_all")

  gt
}

#' @title TrackPlot
#' @description Plot read/UMI coverage and transcript tracks from BAM(s).
#' @param bamfile A path to a BAM file or a list to BAM files. BAM file(s) should be indexed.
#' @param chr Chromosome name.
#' @param start Start position.
#' @param end End position.
#' @param gene Gene name. If provided, chr/start/end are auto-detected from the GTF.
#' @param strand Plot reads on strand. Default is both strands.
#' @param split.bc Split coverage by barcode. Default is FALSE, bulk mode.
#' @param bin Divide plot region into bins. Default is 1000.
#' @param cell.tag Tag for cell barcode in the BAM. Default is "CB".
#' @param umi.tag Tag for UMI in the BAM. Default is "UB".
#' @param gtf GTF database, loaded by gtf2db.
#' @param max.depth Max depth capped to plot. Default is 0 (no capping).
#' @param group.title.size Font size for track group titles.
#' @param cell.group Vector or list of cell group IDs.
#' @param display.genes Vector of genes to display in the gene track.
#' @param meta.features Meta table for feature annotations.
#' @param log.scaled Log-scale the coverage depth. Default FALSE.
#' @param upstream Flank the region upstream. Default 1000.
#' @param downstream Flank the region downstream. Default 1000.
#' @param fragfile Fragment file for ATAC data.
#' @param atac.log.scaled Log-scale ATAC coverage. Default FALSE.
#' @param atac.max.depth Capped depth for ATAC tracks.
#' @param col.by Color BED regions by a metadata column.
#' @param layout.heights Layout heights for the combined plot.
#' @param highlights A region (c(start,end)) or list of regions to highlight.
#' @param junc Also plot junction reads. Default FALSE.
#' @param junc.min.depth Filter junctions below this depth cutoff. Default 5.
#' @importFrom patchwork plot_layout
#' @export
TrackPlot <- function(bamfile=NULL, chr=NULL, start=NULL, end =NULL, gene=NULL,
                      strand = c("both", "forward", "reverse", "ignore"),
                      split.bc = FALSE, bin = 1000, cell.tag = "CB", umi.tag = "UB",
                      gtf = NULL, max.depth = 0, group.title.size = rel(2),
                      cell.group=NULL, display.genes = NULL,  meta.features =NULL,
                      log.scaled = FALSE, upstream = 1000, downstream = 1000,
                      fragfile = NULL,
                      atac.log.scaled = FALSE,
                      atac.max.depth = 0,
                      col.by = NULL, layout.heights =c(1,10,2),
                      highlights = NULL,
                      junc = FALSE, junc.min.depth = 5)
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
  if (end <= start) {
    stop("Corrupt locations (end <= start), check the start and end coordinates.")
  }
  start <- start - upstream
  end <- end + downstream
  message(paste0("Chromosome ", chr, ", start ", start, ", end ", end))
  if (is.null(bamfile)) stop("No bam file specified.")

  strand <- match.arg(strand)
  # ensure layout.heights has at least the 3 elements expected below
  if (length(layout.heights) < 3) {
    layout.heights <- rep(layout.heights, length.out = 3)
  }
  df <- NULL
  if (!is.null(highlights)) {
    if (is.data.frame(highlights)) {
      # data.frame with (xmin, xmax) or (start, end) columns
      if (all(c("xmin", "xmax") %in% colnames(highlights))) {
        df <- highlights[, c("xmin", "xmax")]
      } else if (all(c("start", "end") %in% colnames(highlights))) {
        df <- data.frame(xmin = highlights$start, xmax = highlights$end)
      } else {
        stop("highlights data.frame must have columns 'xmin' and 'xmax', or 'start' and 'end'")
      }
    } else if (is.list(highlights)) {
      # list of vectors: list(c(s1, e1), c(s2, e2))
      highlights <- lapply(highlights, function(x) x[c(1,2)])
      df <- matrix(unlist(highlights), nrow = 2)
      df <- as.data.frame(t(df))
      colnames(df) <- c("xmin", "xmax")
    } else {
      # single vector: c(start, end)
      df <- data.frame(xmin = highlights[1], xmax = highlights[2])
    }

    df <- subset(df, xmin > start & xmax < end)
  }
  
  p1 <- plot.cov(bamfile=bamfile, chr=chr, start=start, end=end, strand=strand, split.bc=split.bc,
                 bin=bin, cell.tag=cell.tag, umi.tag=umi.tag, cell.group=cell.group,
                 log.scaled=log.scaled, max.depth = max.depth, junc=junc, junc.min.depth = junc.min.depth)
  p1 <- p1 + theme(strip.text = element_text(size = group.title.size))
  
  p0 <- NULL
  if (!is.null(meta.features)) {
    if (length(intersect(c("chr","start","end","strand","type"), colnames(meta.features))) != 5)
      stop("No chr/start/end/strand/type column found in meta.features")

    tab <- subset(meta.features, start > 0 & end > 0)
    p0 <- plot.bed(start = start, end = end, peaks = tab, col.by = col.by)
  } 

  p3 <- NULL

  if (!is.null(fragfile)) {
    p3 <- plot.cov2(fragfile = fragfile, chr=chr, start=start, end=end, bin=bin, cell.group=cell.group, log.scaled = atac.log.scaled, max.depth=atac.max.depth)
    p3 <- p3 + theme_cov()
    p3 <- p3 + theme(panel.spacing.y = unit(0.1, "lines"))
    p3 <- p3 + theme(strip.text = element_text(size = group.title.size))
  }
  
  p2 <- plot.genes(chr = chr, start = start, end = end, gtf=gtf, genes=display.genes)
  
  if (!is.null(p0)) {
    if (!is.null(p3)) {
      combined <- p0 / p1 / p3 / p2 + plot_layout(heights=layout.heights[c(1,2,2,3)])
    } else {
      combined <- p0 / p1 / p2 + plot_layout(heights=layout.heights)
    }
  } else if (!is.null(p3)) {
    combined <- p1 / p3 / p2 + plot_layout(heights=layout.heights[c(2,2,3)])
  } else if (length(layout.heights) == 2) {
    combined <- p1 / p2 + plot_layout(heights=layout.heights)
  } else {
    combined <- p1 / p2 + plot_layout(heights=layout.heights[c(2,3)])
  }
  combined <- combined & theme(panel.spacing = unit(0, "lines"), plot.margin = margin(0, 0, 0, 0))

  if (!is.null(df) && nrow(df) > 0L) {
    gt <- patchwork::patchworkGrob(combined)
    gt <- .add_grid_highlights_to_gtable(gt, df, start, end)

    class(gt) <- c("TrackPlot", class(gt))
    gt
  } else {
    combined
  }
}

# S3 print method so TrackPlot objects render correctly
#' @export
#' @method print TrackPlot
print.TrackPlot <- function(x, ...) {
  grid.newpage()
  grid.draw(x)
  invisible(x)
}

# S3 method so grid.draw() works on TrackPlot objects (needed by ggsave)
#' @export
grid.draw.TrackPlot <- function(x, recording = TRUE) {
  class(x) <- setdiff(class(x), "TrackPlot")
  grid.draw(x, recording = recording)
}

#' @importFrom grDevices rgb
#' @importFrom patchwork wrap_plots area plot_layout
#' @importFrom cowplot theme_cowplot
#' @importFrom SeuratObject FetchData
#' @importFrom Seurat SingleDimPlot
#' @importFrom rlang is_integerish
#' @importFrom ggplot2 labs scale_x_continuous scale_y_continuous theme element_rect dup_axis guides element_blank element_text margin scale_color_brewer scale_color_gradientn scale_color_manual coord_fixed ggtitle
RatioPlot0 <- function(object = NULL,
                       assay = NULL,
                       bind.assay = NULL,
                       bind.name = NULL,
                       features = NULL,
                       cells = NULL,
                       dims = c(1,2),
                       cols = c("lightgrey", "red"),
                       pt.size = NULL,
                       alpha = 1,
                       order = TRUE,
                       min.cutoff = NA,
                       max.cutoff = NA,
                       mode = c(1,2,3),
                       reduction = NULL,
                       shape.by = NULL,
                       ncol = NULL,
                       split.by = NULL,
                       by.col = TRUE,
                       coord.fixed = FALSE,
                       combine = TRUE,
                       raster = NULL,
                       raster.dpi = c(512,512),
                       legend.title = NULL,
                       figure_plot = FALSE
                       )
{
  mode <- as.integer(match.arg(as.character(mode), choices = c("1", "2", "3")))
  # Set a theme to remove right-hand Y axis lines
  # Also sets right-hand Y axis text label formatting
  no.right <- theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_text(
      face = "bold",
      size = 14,
      margin = margin(r = 7)
    )
  )
  
  # Get the DimReduc to use
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (!is_integerish(x = dims, n = 2L, finite = TRUE) || !all(dims > 0L)) {
    stop("'dims' must be a two-length integer vector")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% Cells(x = object[[reduction]])
  # Get plotting data
  data <- FetchData(
    object = object,
    vars = c(dims, 'ident'),
    cells = cells
  )

  # Check presence of dimensions
  if (!all(dims %in% colnames(x = data))) {
    stop("The dimensions requested were not found")
  }

  old.assay <- DefaultAssay(object)
  # assay
  DefaultAssay(object) <- assay
  features <- intersect(features, rownames(object))
  if (length(features) == 0) {
    stop("No feature found at assay.")
  }
  
  data1 <- FetchData1(object = object, cells = cells, vars = features, layer = "counts")
  data2 <- NULL
  # binding assay
  if (!is.null(bind.name)) {
    df <- object[[assay]][[]]
    if (bind.name %ni% colnames(df)) {
      stop("No bind.name found in the meta table. Make sure you set the right assay and bind.name.")
    }

    df <- df[features, ]
    features0 <- unique(df[[bind.name]])
    DefaultAssay(object) <- bind.assay
    features0 <- intersect(features0, rownames(object))
    if (length(features0) == 0) {
      stop("No feature found at binding assay.")
    }
    
    data2 <- FetchData1(object = object, cells = cells, vars = features0, layer = "counts")
    data2 <- as.matrix(data2)
    rownames(data2) <- cells
    colnames(data2) <- features0
    df <- df[which(df[[bind.name]] %in% features0),]
    data2 <- as.matrix(data2[, df[[bind.name]]])
    colnames(data2) <- rownames(df)
  } else {
    DefaultAssay(object) <- bind.assay
    features <- intersect(features, rownames(object))
    if (length(features) == 0) {
      stop("No feature found at binding assay.")
    }
    data2 <- FetchData1(object = object, cells = cells, vars = features, layer = "counts")
  }
  # Combine data
  features <- intersect(colnames(data1), colnames(data2))
  if (length(features) == 0) stop("No common features found between the two assays.")

  data1 <- as.matrix(data1[,features])
  data2 <- as.matrix(data2[,features])

  if (mode == 2) {
    data2 <- data2-data1
  } else if (mode == 3) {
    data2 <- data2+data1
  }
  data0 <- data1/data2
  data0[is.infinite(data0)] <- NA
  colnames(data0) <- features
  data <- cbind(data, as.data.frame(data0))
  DefaultAssay(object) <- old.assay

  min.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  )
  max.cutoff <- mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  )

  min.cutoff <- min(min.cutoff)
  max.cutoff <- max(max.cutoff)
  # Figure out splits (FeatureHeatmap)
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(
      EXPR = split.by,
      ident = Idents(object = object)[cells, drop = TRUE],
      object[[split.by, drop = TRUE]][cells, drop = TRUE]
    )
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  
  # Set shaping variable
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  
  # Make list of plots
  plots <- vector(
    mode = "list",
    length = length(x = features) * length(x = levels(x = data$split))  
  )
  # Apply common limits
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))

  # Make the plots
  for (i in 1:length(x = levels(x = data$split))) {
    # Figure out which split we're working with
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    # Make per-feature plots
    for (j in 1:length(x = features)) {
      feature <- features[j]
      cols.use <- NULL
      data.single <- data.plot[, c(dims, 'ident', feature, shape.by)]
      # Make the plot
      plot <- SingleDimPlot(
        data = data.single,
        dims = dims,
        col.by = feature,
        order = order,
        pt.size = pt.size,
        alpha = alpha,
        cols = cols.use,
        shape.by = shape.by,
        label = FALSE,
        raster = raster,
        raster.dpi = raster.dpi
      ) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot() +
        CenterTitle()

      # Make FeatureHeatmaps look nice(ish)
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, colour = 'black'))
        # Add title
        plot <- plot + if (i == 1) {
          if (is.null(legend.title)) {
            labs(title = feature)
          } else {
            labs(title = feature, color = legend.title)
          }
        } else {
          if (is.null(legend.title)) {
            labs(title = feature)
          } else {
            labs(title = NULL, color = legend.title)
          }
        }
        # Add second axis
        if (j == length(x = features)) {
          suppressMessages(
            expr = plot <- plot +
              scale_y_continuous(
                sec.axis = dup_axis(name = ident),
                limits = ylims
              ) +
              no.right
          )
        }
        # Remove left Y axis
        if (j != 1) {
          plot <- plot + theme(
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y.left = element_blank()
          )
        }
        # Remove bottom X axis
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
          )
        }
      } else {
        if (is.null(legend.title)) {
          plot <- plot + labs(title = feature)
        } else {
          plot <- plot + labs(title = feature, color = legend.title)
        }
      }
      # Add colors scale for normal FeaturePlots
      plot <- plot + guides(color = NULL)
      cols.grad <- cols
      if (length(x = cols) > 1) {
        unique.feature.exp <- unique(data.plot[, feature])
        if (length(unique.feature.exp) == 1) {
          warning(paste0(
            "All cells have the same value (",
            unique.feature.exp,
            ") of ",
            dQuote(x = feature)
          ))
          if (unique.feature.exp == 0) {
            cols.grad <- cols[1]
          } else{
            cols.grad <- cols
          }
        }
        plot <- suppressMessages(
          expr = plot + scale_color_gradientn(
            colors = cols.grad,
            guide = "colorbar"
          )
        )
      }
      # Add coord_fixed
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      # Place the plot
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  
  # Remove NULL plots
  plots <- Filter(f = Negate(f = is.null), x = plots)
  # Combine the plots
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) {
      ncol <- 1
    }
    if (length(x = features) > 6) {
      ncol <- 3
    }
    if (length(x = features) > 9) {
      ncol <- 4
    }
  }
  ncol <- ifelse(
    test = is.null(x = split.by),
    yes = ncol, 
    no = length(x = features)
  )
  legend <- split.by %iff% 'none'

  # Transpose the FeatureHeatmap matrix (not applicable for blended FeaturePlots)
  if (isTRUE(x = combine)) {
    if (by.col && !is.null(x = split.by)) {
      plots <- lapply(
        X = plots,
        FUN = function(x) {
          return(suppressMessages(
            expr = x +
              theme_cowplot() +
              ggtitle("") +
              scale_y_continuous(sec.axis = dup_axis(name = ""), limits = ylims) +
              no.right
          ))
        }
      )
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(
          expr = plots[[i]] +
            scale_y_continuous(
              sec.axis = dup_axis(name = features[[idx]]),
              limits = ylims
            ) +
            no.right
        )
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] +
          ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] +
            ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(
        what = rbind,
        args = split(
          x = 1:length(x = plots),
          f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features))
        )
      ))]
      # Set ncol to number of splits (nrow) and nrow to number of features (ncol)
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == 'none') {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol, nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == 'none') {
      plots <- plots & NoLegend()
    }
    plots <- suppressMessages(plots & scale_color_gradientn(colors = cols, limits = c(min.cutoff,max.cutoff)))
  }

  if (isTRUE(figure_plot)) {
    # this parameter is edited from scCustomize::Figure_Plot, credit to original authors
    # See 'samuel-marsh.github.io/scCustomize/articles/FAQ.html' for citation info.
    plots <- plots & NoAxes()
    axis_plot <- ggplot() + geom_segment(data = data.frame(x=c(0,1),y = c(1,0), dx=c(10,1),dy=c(1,10)), aes(x=x,y=y,xend=dx,yend=dy), arrow = grid::arrow(15, grid::unit(0.8, "lines"), ends = "last", type = "closed"))+theme(plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent"), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab(dims[1]) + ylab(dims[2])
    figure_layout <- c(
      area(t = 1, l = 2, b = 11, r = 11),
      area(t = 10, l = 1, b = 12, r = 2))
    plots <- plots + axis_plot + plot_layout(design = figure_layout)
  }
  return(plots)
}

#' This function is edited from Seurat::FeaturePlot, used to visulize the ratio or PSI score on reduction map
#' @title RatioPlot
#' @description Plot ratio score of test feature and its binding feature on reduction map.
#' @param object Seurat object.
#' @param assay Test assay name.
#' @param bind.assay Binding assay name.
#' @param features Features to plot.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character or numeric value corresponding to a palette as specified by RColorBrewer::brewer.pal.info. The default cols is c('lightgrey', 'red'). Set cols = c('lightgrey', 'blue') to get the Seurat 'classical' colors.
#' @param pt.size Adjust point size for plotting
#' @param alpha Alpha value for points
#' @param order Boolean determing whether to plot cells in order of PSI score.
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature
#' @param mode Test mode. For mode 1, X (test feature) vs Y (binding feature). For mode 2, X vs (Y-X). For mode 3, X vs (Y+X). Please note, when set to mode 2 or 3, will use raw counts to update expression value of binding features. Then normalise the counts before testing. For mode 1, will use Layer 'data'. Default is mode 1.
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by A factor in object metadata to split the plot by, pass 'ident' to split by cell identity
#' @param shape.by If NULL, all points are circles (default). You can specify any cell attribute (that can be pulled with FetchData) allowing for both different colors and different shapes on cells.  Only applicable if \code{raster = FALSE}.
#' @param ncol Number of columns to combine multiple features plots to
#' @param coord.fixed Plot cartesian coordinates with fixed aspect ratio
#' @param by.col If splitting by a factor, plot the splits per column with the features as rows 
#' @param combine Combine plots into a single patchwork ggplot object. If \code{FALSE}, return a list of ggplot objects.
#' @param raster If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE if you are encountering that issue (note that plots may take longer to produce/render).
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore(). Default is c(512, 512).
#' @param legend.title Legend title. Default is 'ratio'.
#' @param figure_plot Whether to remove the axes and plot with legend on left of plot denoting axes labels.  (Default is FALSE). This parameter is modified from scCustumize::Figure_Plot, credit to original authors.
#' @return A patchwork ggplot object of \code{combine = TRUE}; otherwise, a list of ggplot objects
#' @export
#' @concept visualization
RatioPlot <- function(object = NULL,
                      assay = NULL,
                      bind.assay = NULL,
                      bind.name = NULL,
                      features = NULL,
                      cells = NULL,
                      dims = c(1,2),
                      cols = c('lightgrey', 'red'),
                      pt.size = NULL,
                      alpha = 1,
                      order = TRUE,
                      min.cutoff = NA,
                      max.cutoff = NA,
                      mode = c(1,2,3),
                      reduction = NULL,
                      shape.by = NULL,
                      ncol = NULL,
                      split.by = NULL,
                      by.col = TRUE,
                      coord.fixed = FALSE,
                      combine = TRUE,
                      raster = NULL,
                      raster.dpi = c(512,512),
                      legend.title = "Ratio",
                      figure_plot = FALSE
                      )
{
  mode <- as.integer(match.arg(as.character(mode), choices = c("1", "2", "3")))
  assay <- assay %||% DefaultAssay(object)
  if (assay %ni% Assays(object)) {
    stop(paste0("Assay ", assay, " is not exist."))
  }

  if (is.null(bind.assay)) {
    stop("bind.assay is not set.")
  }
  if (bind.assay %ni% Assays(object)) {
    stop("bind.assay is not exist, check the binding assay.")
  }

  if (is.null(bind.name)) {
    stop("bind.name is not set.")
  }
  
  plots <- RatioPlot0(object = object,
                      assay = assay,
                      bind.assay = bind.assay,
                      bind.name = bind.name,
                      features = features,
                      cells = cells,
                      dims = dims,
                      cols = cols,
                      pt.size = pt.size,
                      alpha = alpha,
                      order = order,
                      min.cutoff = min.cutoff,
                      max.cutoff = max.cutoff,
                      mode = mode,
                      reduction = reduction,
                      shape.by = shape.by,
                      ncol = ncol,
                      split.by = split.by,
                      by.col = by.col,
                      coord.fixed = coord.fixed,
                      combine = combine,
                      raster = raster,
                      raster.dpi = raster.dpi,
                      legend.title = legend.title,
                      figure_plot = figure_plot)
  
  return(plots)
  
}
#' @title PSIPlot
#' @description Plot PSI score on reduction map.
#' @param object Seurat object.
#' @param exon.assay Exon assay name.
#' @param exclude.assay Excluded exon assay name.
#' @param features Features to plot.
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param cells Vector of cells to plot (default is all cells)
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character or numeric value corresponding to a palette as specified by RColorBrewer::brewer.pal.info. 
#' @param pt.size Adjust point size for plotting
#' @param alpha Alpha value for points
#' @param order Boolean determing whether to plot cells in order of PSI score.
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature
#' @param mode Test mode. For mode 1, X (test feature) vs Y (binding feature). For mode 2, X vs (Y-X). For mode 3, X vs (Y+X). Please note, when set to mode 2 or 3, will use raw counts to update expression value of binding features. Then normalise the counts before testing. For mode 1, will use Layer 'data'. Default is mode 1.
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param group.by Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by A factor in object metadata to split the plot by, pass 'ident' to split by cell identity
#' @param shape.by If NULL, all points are circles (default). You can specify any cell attribute (that can be pulled with FetchData) allowing for both different colors and different shapes on cells.  Only applicable if \code{raster = FALSE}.
#' @param ncol Number of columns to combine multiple features plots to
#' @param coord.fixed Plot cartesian coordinates with fixed aspect ratio
#' @param by.col If splitting by a factor, plot the splits per column with the features as rows 
#' @param combine Combine plots into a single patchwork ggplot object. If \code{FALSE}, return a list of ggplot objects.
# @param raster If true, plot with geom_raster, else use geom_tile. geom_raster may look blurry on some viewing applications such as Preview due to how the raster is interpolated. Set this to FALSE if you are encountering that issue (note that plots may take longer to produce/render).
#' @param raster.dpi Pixel resolution for rasterized plots, passed to geom_scattermore(). Default is c(512, 512).
#' @return A patchwork ggplot object of \code{combine = TRUE}; otherwise, a list of ggplot objects
#' @export
PSIPlot <- function(object = NULL,
                    exon.assay = NULL,
                    exclude.assay = "exclude",
                    features = NULL,
                    dims = c(1,2),
                    cells = NULL,
                    cols = c('lightgrey', 'red'),
                    pt.size = NULL,
                    alpha = 1,
                    order = TRUE,
                    reduction = NULL,
                    shape.by = NULL,
                    ncol = NULL,
                    split.by = NULL,
                    by.col = TRUE,
                    coord.fixed = FALSE,
                    combine = TRUE,
                    raster = NULL,
                    raster.dpi = c(512,512)                    
                    )
{
  exon.assay <- exon.assay %||% DefaultAssay(object)
  if (exon.assay %ni% Assays(object)) {
    stop(paste0("Exon assay ", exon.assay, " is not exist."))
  }

  if (is.null(exclude.assay)) {
    stop("exclude.assay is not set.")
  }
  if (exclude.assay %ni% Assays(object)) {
    stop("exclude.assay is not exist, check the exclude assay.")
  }

  plots <- RatioPlot0(object = object,
                     assay = exon.assay,
                     bind.assay = exclude.assay,
                     features = features,
                     cells = cells,
                     dims = dims,
                     cols = cols,
                     pt.size = pt.size,
                     alpha = alpha,
                     order = order,
                     min.cutoff = 0,
                     max.cutoff = 1,
                     mode = 3,
                     reduction = reduction,
                     shape.by = shape.by,
                     ncol = ncol,
                     split.by = split.by,
                     by.col = by.col,
                     coord.fixed = coord.fixed,
                     combine = combine,
                     raster = raster,
                     raster.dpi = raster.dpi,
                     legend.title = "PSI")
  return(plots)
}
