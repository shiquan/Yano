#' @importFrom reshape2 melt
#' @useDynLib Yano depth2matrix
bamcov0 <- function(bamfile = NULL, chr = NULL, start = -1, end = -1, strand = "both", split.bc = FALSE, cell.group = NULL, bin=1000, cell.tag = "CB", umi.tag = "UB")
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
  
  message(paste0("Process read coverage of ", bamfile))

  strand.flag <- -1
  if (strand == "forward") strand.flag <- 0
  else if(strand == "reverse") strand.flag <- 1
  else if (strand == "ignore") strand.flag <- -2

  if (!is.null(cell.group)) {
    cell.names <- names(cell.group)
    if (is.null(cell.names)) stop("Names of cell.group should not be empty.")
    cell.group <- as.character(cell.group)
    groups <- levels(cell.group) %||% unique(cell.group)
    group.ids <- match(cell.group, groups)
    n <- length(cell.group)
    split.bc <- TRUE    
    dlst <- .Call("depth2matrix", bamfile, chr, start, end, strand.flag, split.bc, cell.tag, umi.tag, 20, cell.names, n, group.ids, groups, FALSE)
  } else {
    dlst <- .Call("depth2matrix", bamfile, chr, start, end, strand.flag, split.bc, cell.tag, umi.tag, 20, NULL, 0, NULL, NULL, FALSE)
  }

  if (is.null(dlst)) return(NULL) ##stop("No cell found in the bam. Check cell names.")
  
  idx.0 <- which(dlst[[3]] == 0)
  idx.1 <- which(dlst[[3]] == 1)

  dlst[[1]] <- as.integer(dlst[[1]]/win)*win

  x <- seq(as.integer(start/win)*win, as.integer(end/win)*win, win)
  y <- unique(dlst[[5]])

  tab <- data.frame()
  
  nr <- length(x)
  nc <- length(y)
  
  if (length(idx.0) > 0) {
    m0 <- Matrix::sparseMatrix(
      i = match(dlst[[1]][idx.0], x),
      j = match(dlst[[5]][idx.0], y),
      x = dlst[[4]][idx.0],
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
      j = match(dlst[[5]][idx.1], y),
      x = dlst[[4]][idx.1],
      dims=c(nr,nc))

    rownames(m1) <- x
    colnames(m1) <- y
    m1 <- as.matrix(m1)
    tab2 <- reshape2::melt(m1)
    tab2$strand <- '-'
    tab <- rbind(tab,tab2)
  }

  colnames(tab) <- c("pos","label","depth","strand")
  tab$depth <- tab$depth/win

  return(tab)
}
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

  if (is.list(bamfile)) {
    nm <- names(bamfile)
    if (is.null(nm)) nm <- c(1:length(bamfile))
    if (is.list(cell.group)) {
      dl <- lapply(nm, function(x) {
        bamcov0(bamfile=bamfile[[x]], chr=as.character(chr), start=start, end=end, strand=strand,
                split.bc=split.bc,
                cell.group=cell.group[[x]], bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)
      })
    } else {
      dl <- lapply(nm, function(x) {
        bamcov0(bamfile=bamfile[[x]], chr=as.character(chr), start=start, end=end, strand=strand,
               split.bc=split.bc,
               cell.group=cell.group, bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)
      })
    }
    
    bc <- bind_rows(dl) %>% group_by(pos, label, strand) %>% summarise(sum(depth, na.rm = TRUE))
    ss <- table(unlist(bc))
    colnames(bc) <- c("pos", "label", "strand", "depth")
  } else {
    bc <- bamcov0(bamfile=bamfile, chr=as.character(chr), start=start, end=end, strand=strand,
                 split.bc=split.bc,
                 cell.group=cell.group, bin=bin, cell.tag=cell.tag, umi.tag=umi.tag)
  }
  
  bc$label <- as.character(bc$label)
  bc$label <- factor(bc$label, levels=gtools::mixedsort(unique(bc$label)))

  return(bc)
}

#'@useDynLib Yano depth2matrix
bamjunc0 <- function(bamfile = NULL, chr = NULL, start = -1, end = -1, strand = "both", split.bc = FALSE, cell.group = NULL, cell.tag = "CB", umi.tag = "UB")
{
  if (is.null(bamfile)) stop("No BAM.")
  if (is.null(chr)) stop("No chromosome name")
  if (start == -1) stop("Start position.")
  if (end == -1) stop("End position")

  len <- end - start
  if (len > 10000000) { ## 10M
    stop("Too large region, try to shrink it before querying.")
  }

  message(paste0("Process junction coverage of ", bamfile))

  strand.flag <- -1
  if (strand == "forward") strand.flag <- 0
  else if(strand == "reverse") strand.flag <- 1
  else if (strand == "ignore") strand.flag <- -2

  if (!is.null(cell.group)) {
    cell.names <- names(cell.group)
    if (is.null(cell.names)) stop("Names of cell.group should not be empty.")
    cell.group <- as.character(cell.group)
    groups <- levels(cell.group) %||% unique(cell.group)
    group.ids <- match(cell.group, groups)
    n <- length(cell.group)
    split.bc <- TRUE    
    dlst <- .Call("depth2matrix", bamfile, chr, start, end, strand.flag, split.bc, cell.tag, umi.tag, 20, cell.names, n, group.ids, groups, TRUE)
  } else {
    dlst <- .Call("depth2matrix", bamfile, chr, start, end, strand.flag, split.bc, cell.tag, umi.tag, 20, NULL, 0, NULL, NULL, TRUE)
  }

  if (is.null(dlst)) return(NULL) ##stop("No cell found in the bam. Check cell names.")
  
  strands <- c("+","-")
  
  tab <- data.frame("start" = dlst[[1]], "end" = dlst[[2]], "label" = dlst[[5]],
                    "depth" = dlst[[4]], "strand" = strands[dlst[[3]]+1])
  
  return(tab)
}
bamjunc <- function(bamfile = NULL, chr = NULL, start = -1, end = -1, strand = "both", split.bc = FALSE, cell.group = NULL, cell.tag = "CB", umi.tag = "UB")
{
  if (is.null(bamfile)) stop("No BAM.")
  if (is.null(chr)) stop("No chromosome name")
  if (start == -1) stop("Start position.")
  if (end == -1) stop("End position")

  len <- end - start
  if (len > 10000000) { ## 10M
    stop("Too large region, try to shrink it before querying.")
  }

  if (is.list(bamfile)) {
    nm <- names(bamfile)
    if (is.list(cell.group)) {
      dl <- lapply(nm, function(x) {
        bamjunc0(bamfile=bamfile[[x]], chr=as.character(chr), start=start, end=end, strand=strand,
                 split.bc=split.bc,
                 cell.group=cell.group[[x]], cell.tag=cell.tag, umi.tag=umi.tag)
      })
    } else {
      dl <- lapply(nm, function(x) {
        bamjunc0(bamfile=bamfile[[x]], chr=as.character(chr), start=start, end=end, strand=strand,
                 split.bc=split.bc,
                 cell.group=cell.group, cell.tag=cell.tag, umi.tag=umi.tag)
      })
    }
    
    juncs <- bind_rows(dl) %>%  group_by(start, end , label, strand) %>% summarise(sum(depth, na.rm = TRUE))
    colnames(juncs) <- c("start", "end", "label", "strand", "depth")

  } else {
    juncs <- bamjunc0(bamfile=bamfile, chr=as.character(chr), start=start, end=end, strand=strand,
                      split.bc=split.bc,
                      cell.group=cell.group, cell.tag=cell.tag, umi.tag=umi.tag)
  }
  juncs$label <- as.character(juncs$label)
  juncs$label <- factor(juncs$label, levels=gtools::mixedsort(unique(juncs$label)))

  return(juncs)
}
fragcov0 <- function(fragfile = NULL, chr = NULL, start = -1, end = -1, split.bc = FALSE, cell.group = NULL, bin=1000)
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

  if (is.list(fragfile)) {
    nm <- names(fragfile)
    if (is.list(cell.group)) {
      dl <- lapply(nm, function(x) {
        fragcov0(fragfile=fragfile[[x]], chr=as.character(chr), start=start, end=end,
                split.bc=split.bc,cell.group=cell.group[[x]], bin=bin)
      })
    } else {
      dl <- lapply(nm, function(x) {
        fragcov0(fragfile=fragfile[[x]], chr=as.character(chr), start=start, end=end,
                split.bc=split.bc, cell.group=cell.group, bin=bin)
      })
    }

    bc <- bind_rows(dl) %>%  group_by(pos, label, strand) %>% summarise(sum(depth, na.rm = TRUE))
    ss <- table(unlist(bc))
    colnames(bc) <- c("pos", "label", "strand", "depth")
  } else {
    bc <- fragcov0(fragfile=fragfile, chr=as.character(chr), start=start, end=end,
                  split.bc=split.bc, cell.group=cell.group, bin=bin)
  }
  
  bc$label <- as.character(bc$label)
  bc$label <- factor(bc$label, levels=gtools::mixedsort(unique(bc$label)))
  return(bc)
}
