# =========================================================================
# Internal helpers -- shared across all selector functions
# =========================================================================

# Robust cell extraction: get selectable cell names from a ggplot object.
# Falls back to colnames(object) if rownames are absent, warns on partial
# mismatches instead of silently selecting all cells.
.extract_plot_cells <- function(p, object) {
  cells <- rownames(p$data)
  if (is.null(cells) || length(cells) == 0) {
    return(colnames(object))
  }
  bad <- setdiff(cells, colnames(object))
  if (length(bad) > 0) {
    if (length(bad) == length(cells)) {
      warning("No plot rownames match object column names; using all cells.")
      return(colnames(object))
    }
    warning("Dropping ", length(bad), " cell name(s) not present in object.")
    cells <- intersect(cells, colnames(object))
  }
  if (length(cells) == 0) stop("No selectable cells found in plot data.")
  cells
}

# Internal: concave-hull filter for spatial/dim reduction selections
.concave_filter <- function(p, idx, data_cols = c("x", "y")) {
  data <- p$data[idx, data_cols, drop = FALSE]
  data[] <- lapply(data, as.numeric)
  data <- as.matrix(data)
  keys <- paste0(round(data[, 1], 4), "_", round(data[, 2], 4))
  ch <- concaveman::concaveman(data)
  ch_keys <- paste0(round(ch[, 1], 4), "_", round(ch[, 2], 4))
  idx[keys %in% ch_keys]
}

# Internal: mark selected cells and re-plot
.plot_selected <- function(object, cells, plot_func, ...) {
  object[["selector"]] <- "unselect"
  object[["selector"]][cells, ] <- "selected"
  p <- plot_func(object, group.by = "selector", ...)
  print(p)
}

# =========================================================================
# Core engine -- shared by all selector functions
# =========================================================================

# Build a plot, let the user interactively select cells, optionally apply
# a concave-hull filter, and return either the cell names or a subset object.
.select_cells <- function(object,
                          plot_func,
                          ...,
                          return.object = FALSE,
                          plot.selected = TRUE,
                          concave.hull = FALSE,
                          p = NULL) {
  # 1. Build the plot if not supplied, then draw it
  if (is.null(p)) {
    p <- plot_func(object, ..., combine = FALSE)
    if (is.list(p)) p <- p[[1L]]
  }
  grid.draw(p)

  # 2. Get selectable cell names
  cells <- .extract_plot_cells(p, object)

  # 3. Interactive selection via the ggplot graphics device
  rl <- ggplot_selector()
  idx <- rl[["idx"]]
  cells <- cells[idx]

  # 4. Optional concave-hull refinement
  if (concave.hull && length(cells) > 0) {
    keep <- .concave_filter(p, idx)
    cells <- cells[keep]
    cells <- as.character(cells)
  }

  # 5. Optional re-plot with selection highlighted
  if (plot.selected) {
    .plot_selected(object, cells, plot_func, ...)
  }

  if (return.object) object[, cells] else cells
}

# =========================================================================
# Public selector functions -- thin wrappers around .select_cells
# =========================================================================

#' @title DimSelector
#' @description Select cells from dimension plot.
#' @param object Seurat object.
#' @param return.object Return Seurat object. Default returns selected cell IDs.
#' @param plot.selected Plot the selected cells.
#' @param ... Parameters passed to Seurat::DimPlot().
#' @importFrom Seurat DimPlot
#' @importFrom grid grid.draw
#' @export
DimSelector <- function(object,
                        return.object = FALSE,
                        plot.selected = TRUE,
                        combine = FALSE,
                        ...) {
  .select_cells(object,
                plot_func = DimPlot,
                return.object = return.object,
                plot.selected = plot.selected,
                combine = combine,
                ...)
}

#' @title FeatureSelector
#' @description Select cells from Feature plot.
#' @param object Seurat object.
#' @param feature Feature to plot.
#' @inheritParams DimSelector
#' @export
FeatureSelector <- function(object,
                            feature,
                            return.object = FALSE,
                            plot.selected = TRUE,
                            combine = FALSE,
                            reduction = NULL,
                            cells = NULL,
                            ...) {
  if (length(feature) == 0) stop("No feature specified.")
  if (length(feature) > 1) {
    warning("Only one feature is supported. Using the first feature.")
    feature <- feature[1L]
  }
  # pre-build the FeaturePlot to extract data for cell-name fallback
  p <- FeaturePlot(object, features = feature, ..., combine = FALSE)
  if (is.list(p)) p <- p[[1L]]
  .select_cells(object,
                plot_func = FeaturePlot,
                p = p,
                features = feature,
                return.object = return.object,
                plot.selected = plot.selected,
                combine = combine,
                reduction = reduction,
                cells = cells,
                ...)
}

#' @title SpatialSelector
#' @description Select cells from Spatial plot.
#' @inheritParams DimSelector
#' @param group.by Grouping variable.
#' @param ... Parameters passed to Seurat::SpatialPlot().
#' @importFrom Seurat SpatialPlot
#' @importFrom grid grid.draw
#' @export
SpatialSelector <- function(object,
                            group.by = NULL,
                            return.object = FALSE,
                            plot.selected = TRUE,
                            ...) {
  .select_cells(object,
                plot_func = SpatialPlot,
                group.by = group.by,
                return.object = return.object,
                plot.selected = plot.selected,
                ...)
}

#' @title ImageDimSelector
#' @description Select cells from Image dimension plot.
#' @inheritParams DimSelector
#' @param group.by Grouping variable.
#' @param ... Parameters passed to Seurat::ImageDimPlot().
#' @importFrom Seurat ImageDimPlot
#' @importFrom grid grid.draw
#' @export
ImageDimSelector <- function(object,
                             group.by = NULL,
                             return.object = FALSE,
                             plot.selected = TRUE,
                             combine = FALSE,
                             ...) {
  .select_cells(object,
                plot_func = ImageDimPlot,
                group.by = group.by,
                return.object = return.object,
                plot.selected = plot.selected,
                combine = combine,
                ...)
}

#' @title SpatialConcaveHull
#' @description Select cells from spatial plot using concave hull lasso.
#' @inheritParams DimSelector
#' @param group.by Grouping variable.
#' @param knn Number of nearest neighbors for concave hull. Default is 3.
#' @param ... Parameters passed to Seurat::SpatialPlot().
#' @importFrom Seurat SpatialPlot
#' @importFrom grid grid.draw
#' @export
SpatialConcaveHull <- function(object,
                               group.by = NULL,
                               return.object = FALSE,
                               plot.selected = TRUE,
                               knn = 3,
                               ...) {
  .select_cells(object,
                plot_func = SpatialPlot,
                group.by = group.by,
                return.object = return.object,
                plot.selected = plot.selected,
                concave.hull = TRUE,
                ...)
}

#' @title DimConcaveHull
#' @description Select cells from dimension reduction plot using concave hull lasso.
#' @inheritParams DimSelector
#' @param ... Parameters passed to Seurat::DimPlot().
#' @importFrom Seurat DimPlot
#' @importFrom grid grid.draw
#' @export
DimConcaveHull <- function(object,
                           return.object = FALSE,
                           plot.selected = TRUE,
                           ...) {
  .select_cells(object,
                plot_func = DimPlot,
                return.object = return.object,
                plot.selected = plot.selected,
                concave.hull = TRUE,
                ...)
}

# =========================================================================
# Interactive selection engine -- device-level lasso tool
# =========================================================================

#' @importFrom ggplot2 ggplot_build ggplot_gtable
#' @importFrom stringr str_detect
#' @importFrom grid deviceLoc seekViewport
#' @title ggplot_selector
#' @description Interactive cell selection tool for ggplot2-based plots.
#' @return A list with x and y coordinates and logical index of selected cells.
#' @export
ggplot_selector <- function() {
  if (.Platform$GUI == "RStudio" &&
      str_detect(.Device, "(RStudio)|(null device)")) {
    stop("Selector does not work in the RStudio plot viewer. ",
         "Create a different device with grDevices::x11 (Linux) ",
         "or grDevices::quartz (Mac OS).")
  }

  plot_info <- ggplot_build(last_plot())
  gtable <- ggplot_gtable(plot_info)

  panel_info <- gtable$layout[grepl("^panel", gtable$layout$name), ]

  if (nrow(panel_info) == 0) {
    stop("No panel found. Print your plot first.")
  }

  panel_info <- panel_info[, c("t", "l", "b", "r", "name")]
  panel_info$vp_name <- with(panel_info,
    sprintf("%s.%d-%d-%d-%d", name, t, r, b, l))

  seekViewport(panel_info$vp_name)
  on.exit(try(upViewport(), silent = TRUE), add = TRUE)

  xrange <- plot_info$layout$panel_params[[1]]$x.range
  yrange <- plot_info$layout$panel_params[[1]]$y.range

  data <- plot_info$data[[1]]

  x0 <- (data$x - xrange[1]) / diff(xrange)
  y0 <- (data$y - yrange[1]) / diff(yrange)

  sl <- deviceLoc(unit(x0, "npc"), unit(y0, "npc"))

  ta <- .Call("selector")
  n <- ta[[3]]
  if (n == 0) {
    stop("No location selected.")
  }
  if (n < 3) {
    stop("Less than 3 edges, cannot select any point.")
  }

  rl <- list()
  seekViewport(panel_info$vp_name)
  x1 <- as.numeric(convertX(unit(ta[[1]], "inches"), unitTo = "npc"))
  y1 <- as.numeric(convertY(unit(ta[[2]], "inches"), unitTo = "npc"))
  rl[["x"]] <- x1 * diff(xrange) + xrange[1]
  rl[["y"]] <- y1 * diff(yrange) + yrange[1]

  idx <- .Call("points_in_polygen", ta[[1]], ta[[2]],
               as.numeric(sl$x), as.numeric(sl$y))
  rl[["idx"]] <- idx
  rl
}
