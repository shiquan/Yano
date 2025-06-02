#' @title DimSelector
#' @description Select cells from dimension plot. Parameters most inhert from Seurat::DimPlot().
#' 
#' @importFrom Seurat DimPlot
#' @export
DimSelector <- function(object,
                        dims = c(1, 2),
                        cells = NULL,
                        cols = NULL,
                        pt.size = NULL,
                        reduction = NULL,
                        group.by = NULL,
                        shape.by = NULL,
                        order = NULL,
                        shuffle = FALSE,
                        seed = 1,
                        label = FALSE,
                        label.size = 4,
                        label.color = 'black',
                        label.box = FALSE,
                        repel = FALSE,
                        alpha = 1,
                        cells.highlight = NULL,
                        cols.highlight = '#DE2D26',
                        sizes.highlight = 1,
                        na.value = 'grey50',
                        raster = NULL,
                        raster.dpi = c(512, 512),
                        return.object = FALSE,
                        plot.selected = TRUE
                        )
{
  p <- DimPlot(object,
               dims = dims,
               cells = cells,
               cols = cols,
               pt.size = pt.size,
               reduction = reduction,
               group.by = group.by,
               shape.by = shape.by,
               order = order,
               shuffle = shuffle,
               seed = seed,
               label = label,
               label.size = label.size,
               label.color = label.color,
               label.box = label.box,
               repel = repel,
               alpha = alpha,
               cells.highlight = cells.highlight,
               cols.highlight = cols.highlight,
               sizes.highlight = sizes.highlight,
               na.value = na.value,
               raster = raster,
               raster.dpi = raster.dpi,
               combine=FALSE
               )
  if (is.list(p)) {
    p <- p[[1]]
  }

  grid.draw(p)

  cells0 <- rownames(p$data)
  if (length(intersect(cells0, colnames(object))) != length(cells0)) {
    cells0 <- colnames(object)
  }
    
  idx <- ggplot_selector()
  
  cells0 <- cells0[which(idx)]
  
  if (plot.selected) {
    object[['selector']] <- 'unselect'
    object[['selector']][cells0,] <- 'selected'
    
    p <- DimPlot(object, group.by = 'selector',
                 dims = dims,
                 cells = cells,
                 cols = cols,
                 pt.size = pt.size,
                 reduction = reduction,
                 shape.by = shape.by,
                 order = order,
                 shuffle = shuffle,
                 seed = seed,
                 label = label,
                 label.size = label.size,
                 label.color = label.color,
                 label.box = label.box,
                 repel = repel,
                 alpha = alpha,
                 cells.highlight = cells.highlight,
                 cols.highlight = cols.highlight,
                 sizes.highlight = sizes.highlight,
                 na.value = na.value,
                 raster = raster,
                 raster.dpi = raster.dpi)
    print(p)
  }

  if (return.object) {
    object[,cells0]
  } else {
    cells0
  }  
}

#' @export
FeatureSelector <- function(object,
                            feature,
                            dims = c(1, 2),
                            cells = NULL,
                            cols = if (blend) {
                              c('lightgrey', '#ff0000', '#00ff00')
                            } else {
                              c('lightgrey', 'blue')
                            },
                            pt.size = NULL,
                            alpha = 1,
                            order = FALSE,
                            min.cutoff = NA,
                            max.cutoff = NA,
                            reduction = NULL,
                            keep.scale = "feature",
                            shape.by = NULL,
                            slot = 'data',
                            blend = FALSE,
                            blend.threshold = 0.5,
                            label = FALSE,
                            label.size = 4,
                            label.color = "black",
                            repel = FALSE,
                            ncol = NULL,
                            coord.fixed = FALSE,
                            raster = NULL,
                            raster.dpi = c(512, 512),
                            return.object = FALSE,
                            plot.selected = TRUE

                            )
{
  if (length(feature) == 0) {
    stop("No feature specified.")
  }
  if (length(feature) > 1) {
    warning("Only one feature is support. Here only use the first feature.")
    feature <- feature[1]
  }
  p <- FeaturePlot(object,
                   features = feature,
                   dims = dims,
                   cells = cells,
                   cols = cols,
                   pt.size = pt.size,
                   alpha = alpha,
                   order = order,
                   min.cutoff = min.cutoff,
                   max.cutoff = max.cutoff,
                   reduction = reduction,
                   keep.scale = keep.scale,
                   shape.by = shape.by,
                   slot = slot,
                   blend = blend,
                   blend.threshold = blend.threshold,
                   label = label,
                   label.size = label.size,
                   label.color = label.color,
                   repel = repel,
                   coord.fixed = coord.fixed,
                   raster = raster,
                   raster.dpi = c(512, 512),
                   combine = FALSE)
  if (is.list(p)) {
    p <- p[[1]]
  }

  grid.draw(p)

  cells0 <- rownames(p$data)
  if (length(intersect(cells0, colnames(object))) != length(cells0)) {
    cells0 <- colnames(object)
  }
    
  idx <- ggplot_selector()
  
  cells0 <- cells0[which(idx)]
  
  if (plot.selected) {
    object[['selector']] <- 'unselect'
    object[['selector']][cells0,] <- 'selected'
    p <- DimPlot(object, group.by = 'selector', reduction = reduction, cells =cells) 
    print(p)
  }

  if (return.object) {
    object[,cells0]
  } else {
    cells0
  }
  
}
#' @importFrom ggplot2 scale_fill_gradientn ggtitle theme element_text ggplot_build last_plot
#' @importFrom Seurat Cells Images DefaultAssay Idents SpatialPlot
#' @importFrom grid current.vpTree current.vpPath seekViewport upViewport current.transform convertX convertY unit grid.draw
#' @importFrom scales hue_pal
#' @export
SpatialSelector <- function(
  object,
  group.by = NULL,
  feature = NULL,
  image = NULL,
  cols = NULL,
  image.alpha = 1,
  image.scale = "lowres",
  crop = TRUE,
  slot = 'data',
  label = FALSE,
  label.size = 5,
  label.color = 'white',
  label.box = TRUE,
  repel = FALSE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  shape = 21,
  return.object = FALSE,
  plot.selected = TRUE
  ) {
  if (length(feature) > 1) {
    warning("SpatialSelector only support plot one feature.")
    feature <- feature[1]
  }
  p <- SpatialPlot(object,
                   group.by = group.by,
                   features = feature,
                   images = image,
                   cols = cols,
                   image.alpha = image.alpha,
                   image.scale = image.scale,
                   crop = crop,
                   slot = slot,
                   label = label,
                   label.size = label.size,
                   label.color = label.color,
                   repel = repel,
                   pt.size.factor = pt.size.factor,
                   alpha = alpha,
                   shape = shape,
                   combine = FALSE
                   )
  if (is.list(p)) {
    p <- p[[1]]
  }
  grid.draw(p)
  

  cells <- rownames(p$data)
  if (length(intersect(cells, colnames(object))) != length(cells)) {
    cells <- colnames(object)
  }
    
  idx <- ggplot_selector()
  
  cells <- cells[which(idx)]
  
  if (plot.selected) {
    object[['selector']] <- 'unselect'
    object[['selector']][cells,] <- 'selected'
    p <- SpatialPlot(object, group.by = 'selector', images = image)
    print(p)
  }

  if (return.object) {
    object[,cells]
  } else {
    cells
  }
}

#' @importFrom Seurat ImageDimPlot
#' @export
ImageDimSelector <- function(object,
                             fov = NULL,
                             boundaries = NULL,
                             group.by = NULL,
                             split.by = NULL,
                             cols = NULL,
                             shuffle.cols = FALSE,
                             size = 0.5,
                             molecules = NULL,
                             mols.size = 0.1,
                             mols.cols = NULL,
                             mols.alpha = 1.0,
                             nmols = 1000,
                             alpha = 1.0,
                             border.color = 'white',
                             border.size = NULL,
                             na.value = 'grey50',
                             dark.background = TRUE,
                             crop = FALSE,
                             cells = NULL,
                             overlap = FALSE,
                             axes = FALSE,
                             combine = TRUE,
                             coord.fixed = TRUE,
                             flip_xy = TRUE,
                             return.object = FALSE,
                             plot.selected = TRUE
                             )
{
  p <- ImageDimPlot(object,
                    fov=fov,
                    boundaries = boundaries,
                    group.by = group.by,
                    cols = cols,
                    shuffle.cols = shuffle.cols,
                    size = size,
                    molecules = molecules,
                    mols.size = mole.size,
                    mols.cols = mols.cols,
                    mols.alpha = mols.alpha,
                    nmols = nmols,
                    alpha = alpha,
                    border.color = border.color,
                    border.size = border.size,
                    na.value = na.value,
                    dark.background = dark.background,
                    crop = crop,
                    cells = cells,
                    overlap = overlap,
                    axes = axes,
                    coord.fixed = coord.fixed,
                    flip_xy = flip_xy,
                    combine = FALSE)

  if (is.list(p)) {
    p <- p[[1]]
  }
  grid.draw(p)
  data <- p$data
  if ("cell" %in% colnames(p$data)) {
    cells0 <- data$cell
    cells0 <- gsub("centroids_", "", cells0)
  } else {
    cells0 <- rownames(data)
  }
  if (length(intersect(cells0, colnames(object))) != length(cells0)) {
    cells0 <- colnames(object)
  }
    
  idx <- ggplot_selector()
  
  cells0 <- cells0[which(idx)]
  
  if (plot.selected) {
    object[['selector']] <- 'unselect'
    object[['selector']][cells0,] <- 'selected'
    
    p <- ImageDimPlot(object, group.by = 'selector',
                      fov=fov,
                      boundaries = boundaries,
                      cols = cols,
                      shuffle.cols = shuffle.cols,
                      size = size,
                      molecules = molecules,
                      mols.size = mole.size,
                      mols.cols = mols.cols,
                      mols.alpha = mols.alpha,
                      nmols = nmols,
                      alpha = alpha,
                      border.color = border.color,
                      border.size = border.size,
                      na.value = na.value,
                      dark.background = dark.background,
                      crop = crop,
                      cells = cells,
                      overlap = overlap,
                      axes = axes,
                      coord.fixed = coord.fixed,
                      flip_xy = flip_xy)
    print(p)
  }

  if (return.object) {
    object[,cells0]
  } else {
    cells0
  }
}

#' @importFrom Seurat SpatialPlot
#' @export
SpatialConcaveHull <- function(
  object,
  group.by = NULL,
  feature = NULL,
  image = NULL,
  cols = NULL,
  image.alpha = 1,
  image.scale = "lowres",
  crop = TRUE,
  slot = 'data',
  label = FALSE,
  label.size = 5,
  label.color = 'white',
  label.box = TRUE,
  repel = FALSE,
  pt.size.factor = 1.6,
  alpha = c(1, 1),
  shape = 21,
  return.object = FALSE,
  plot.selected = TRUE,
  knn = 3
  ) {
  if (length(feature) > 1) {
    warning("SpatialSelector only support plot one feature.")
    feature <- feature[1]
  }
  p <- SpatialPlot(object,
                   group.by = group.by,
                   features = feature,
                   images = image,
                   cols = cols,
                   image.alpha = image.alpha,
                   image.scale = image.scale,
                   crop = crop,
                   slot = slot,
                   label = label,
                   label.size = label.size,
                   label.color = label.color,
                   repel = repel,
                   pt.size.factor = pt.size.factor,
                   alpha = alpha,
                   shape = shape,
                   combine = FALSE
                   )
  if (is.list(p)) {
    p <- p[[1]]
  }
  grid.draw(p)
  

  cells <- rownames(p$data)
  if (length(intersect(cells, colnames(object))) != length(cells)) {
    cells <- colnames(object)
  }
    
  idx <- ggplot_selector()

  cells <- cells[which(idx)]
  data <- p$data[which(idx),c("x", "y")]
  data$x <- as.numeric(data$x)
  data$y <- as.numeric(data$y)
  data <- as.matrix(data)
  cells0 <- paste0(round(data[,1],4),"_",round(data[,2],4))
  names(cells) <- cells0
  ch <- concavemn::concaveman(data)
  ch.ps <- paste0(ch[,1], "_", ch[,2])
  cells <- cells[ch.ps]
  cells <- as.character(cells)

  ## ta <- .Call("concave_hull", data$x, data$y)
  ## cells <- cells[ta[[3]]]
  
  if (plot.selected) {
    object[['selector']] <- 'unselect'
    object[['selector']][cells,] <- 'selected'
    
    p <- SpatialPlot(object, group.by = 'selector', images = image)
    print(p)
  }

  if (return.object) {
    object[,cells]
  } else {
    cells
  }
}

#' @importFrom Seurat DimPlot
#' @export
DimConcaveHull <- function(
                           object,
                           dims = c(1, 2),
                           cells = NULL,
                           cols = NULL,
                           pt.size = NULL,
                           reduction = NULL,
                           group.by = NULL,
                           shape.by = NULL,
                           order = NULL,
                           shuffle = FALSE,
                           seed = 1,
                           label = FALSE,
                           label.size = 4,
                           label.color = 'black',
                           label.box = FALSE,
                           repel = FALSE,
                           alpha = 1,
                           cells.highlight = NULL,
                           cols.highlight = '#DE2D26',
                           sizes.highlight = 1,
                           na.value = 'grey50',
                           raster = NULL,
                           raster.dpi = c(512, 512),
                           return.object = FALSE,
                           plot.selected = TRUE
                           )
{
  p <- DimPlot(object,
               dims = dims,
               cells = cells,
               cols = cols,
               pt.size = pt.size,
               reduction = reduction,
               group.by = group.by,
               shape.by = shape.by,
               order = order,
               shuffle = shuffle,
               seed = seed,
               label = label,
               label.size = label.size,
               label.color = label.color,
               label.box = label.box,
               repel = repel,
               alpha = alpha,
               cells.highlight = cells.highlight,
               cols.highlight = cols.highlight,
               sizes.highlight = sizes.highlight,
               na.value = na.value,
               raster = raster,
               raster.dpi = raster.dpi,
               combine=FALSE
               )
  if (is.list(p)) {
    p <- p[[1]]
  }

  grid.draw(p)

  cells0 <- rownames(p$data)
  if (length(intersect(cells0, colnames(object))) != length(cells0)) {
    cells0 <- colnames(object)
  }
    
  idx <- ggplot_selector()

  cells0 <- cells0[which(idx)]
  data <- p$data[which(idx),c("x", "y")]
  data$x <- as.numeric(data$x)
  data$y <- as.numeric(data$y)
  data <- as.matrix(data)
  cells1 <- paste0(round(data[,1],4),"_",round(data[,2],4))
  names(cells0) <- cells1
  ch <- concaveman::concaveman(data)
  ch.ps <- paste0(ch[,1], "_", ch[,2])
  cells0 <- cells0[ch.ps]
  cells0 <- as.character(cells0)

  
  if (plot.selected) {
    object[['selector']] <- 'unselect'
    object[['selector']][cells0,] <- 'selected'
    
    p <- DimPlot(object, group.by = 'selector', reduction = reduction, cells=cells)
    print(p)
  }

  if (return.object) {
    object[,cells0]
  } else {
    cells0
  }
}

#' @importFrom ggplot2 ggplot_build ggplot_gtable
#' @importFrom stringr str_detect
#' @importFrom grid deviceLoc seekViewport
#' @export
ggplot_selector <- function()
{
  if (
    .Platform$GUI == "RStudio" && str_detect(.Device, "(RStudio)|(null device)")
  ) {
    stop(
      "Now selector not works in the RStudio plot viewer. Create a different device with grDevices::x11 (Linux) or grDevices::quartz (Mac OS)."
    )
  }

  plot_info <- ggplot_build(last_plot())
  gtable <- ggplot_gtable(plot_info)

  panel_info <- gtable$layout[grepl("^panel", gtable$layout$name),]

  if (length(panel_info) == 0) {
    stop("Print your plot first.")
  }
  
  panel_info <- panel_info[, c("t", "l", "b", "r", "name")]
  panel_info$vp_name <- with(panel_info, sprintf("%s.%d-%d-%d-%d", name, t, r, b, l))

  seekViewport(panel_info$vp_name)

  xrange <- plot_info$layout$panel_params[[1]]$x.range
  yrange <- plot_info$layout$panel_params[[1]]$y.range

  data <- plot_info$data[[1]]

  x0 <- (data$x - xrange[1])/diff(xrange)
  y0 <- (data$y - yrange[1])/diff(yrange)

  sl <- deviceLoc(unit(x0,"npc"),unit(y0, "npc"))
  
  ta <- .Call('selector')
  n <- ta[[3]]
  if (n == 0) {
    stop("No select any location.")
  }

  if (n < 3) {
    stop("Less than 3 edges, cannot select any point.")
  }
  
  idx <- .Call("points_in_polygen", ta[[1]], ta[[2]], as.numeric(sl$x), as.numeric(sl$y))
  
  idx
}
