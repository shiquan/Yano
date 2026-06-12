#' @title geom_splice
#' @description A ggplot2 geom for drawing splice junction arcs.
#' @inheritParams ggplot2::layer
#' @param spline_shape Spline shape parameter, default -0.25.
#' @export
geom_splice <- function(mapping = NULL,
                        data = NULL, 
                        stat = "splice", 
                        position = "identity", 
                        ..., 
                        spline_shape = -0.25,
                        lineend = "butt", 
                        na.rm = FALSE, 
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping, 
    stat = stat,
    geom = GeomPath,
    position = position, 
    show.legend = show.legend, 
    inherit.aes = inherit.aes, 
    params = list(
      lineend = lineend, 
      na.rm = na.rm,
      spline_shape = spline_shape,
      ...
    )
  )
}

## #' @importFrom ggplot2 ggproto
## #' @importFrom grid gpar
## #' @export
## GeomSplice <- ggproto("GeomSplice", Geom,
##                       draw_panel = function(data,
##                                             panel_params,
##                                             coord,
##                                             lineend = "butt",
##                                             spline_shape = -0.25,
##                                             na.rm = FALSE) {
                        
##                         data <- remove_missing(
##                           df = data,
##                           na.rm = na.rm,
##                           vars = c("x", "y", "xend", "height"),
##                           name = "geom_splice"
##                         )
                        
##                         if (is.null(data) || nrow(data) == 0) return(zeroGrob())
##                         if (!coord$is_linear()) {
##                           rlang::warn(
##                             "spline geom only works correctly on linear coordinate systems"
##                           )
##                         }
##                         data <- coord$transform(data, panel_params)
    
##                         spliceGrob(
##                           data = data,
##                           #default.units = "npc",
##                           spline_shape = spline_shape,
##                           gp = gpar(
##                             col = alpha(data$colour, data$alpha),
##                             lwd = data$linewidth * .pt,
##                             lty = data$linetype,
##                             lineend = lineend
##                           )
##                         )
##                       },
                      
##                       required_aes = c("x", "y", "xend", "height"),
##                       default_aes = aes(
##                         colour = "black", 
##                         linewidth = 0.5, 
##                         linetype = 1L, 
##                         alpha = NA
##                       )
##                       )


#' @title stat_splice
#' @description A ggplot2 stat for computing splice junction paths.
#' @inheritParams ggplot2::layer
#' @param spline_shape Spline shape parameter, default -0.25.
#' @param open Whether the spline is open, default TRUE.
#' @param rep_ends Whether to repeat ends, default TRUE.
#' @export
stat_splice <- function(mapping = NULL, data = NULL, geom = "splice",
                        position = "identity", na.rm = TRUE, show.legend = NA, inherit.aes = TRUE,
                        spline_shape=-0.25, open=TRUE, rep_ends=TRUE, ...) {
  layer(
    stat = StatSplice,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(spline_shape=spline_shape,
                  open=open,
                  na.rm = na.rm,
                  rep_ends=rep_ends,
                  ...
    )
  )
}

#' @title StatSplice
#' @description A ggplot2 Stat for computing splice junction paths.
#' @export
#' @importFrom graphics xspline
StatSplice <- ggproto("StatSplice", Stat,
                       setup_data = function(data, params) {
                         if (anyDuplicated(data$group)) {
                           data$group <- paste(data$group, seq_len(nrow(data)), sep = "-")
                         }
                         data
                       },
                      
                      compute_panel = function(data, scales, spline_shape=-0.25) {
                        cols_to_keep <- setdiff(names(data), c("x", "y", "xend", "height"))
                        splices <- lapply(seq_len(nrow(data)), function(i) {
                          splice_path <- create_splice(
                            data$x[i], 
                            data$y[i], 
                            data$xend[i], 
                            data$height[i], 
                            spline_shape
                          )
                          cbind(splice_path, unclass(data[i, cols_to_keep]))
                        })
                        do.call(rbind, splices)
                      },
                      required_aes = c("x", "y", "xend", "height")
                      )

#' @title create_splice
#' @description Create a spline path for a splice junction arc.
#' @param x Start x-coordinate.
#' @param y Start y-coordinate.
#' @param xend End x-coordinate.
#' @param height Height of the arc.
#' @param spline_shape Spline shape parameter, default -0.25.
#' @return A data.frame with x and y coordinates for the spline path.
#' @importFrom graphics xspline
#' @export
create_splice <- function(x, y, xend, height = 0, spline_shape = -0.25)
{
  xm <- (x+xend)/2
  ym <- y+height
  # xspline needs an initialized graphics device even with draw=FALSE;
  # use pdf(file=NULL) to avoid writing to disk.  on.exit guarantees
  # the temporary device is cleaned up even if xspline throws an error.
  pdf(file = NULL)
  on.exit(invisible(dev.off()), add = TRUE, after = FALSE)
  plot.new()
  tmp <- xspline(x=c(x,xm,xend), y = c(y,ym,y), spline_shape, TRUE, TRUE, draw = FALSE)
  data.frame(x=tmp$x, y=tmp$y)
}

## #' @importFrom grid unit is.unit gTree
## #' @export
## spliceGrob <- function(data,
##                        #default.units = "native",
##                        spline_shape = -0.25,
##                        name = NULL,
##                        gp = gpar(), 
##                        vp = NULL) {
##   gTree(
##     data = data,
##     spline_shape = spline_shape,
##     name = name,
##     gp = gp,
##     vp = vp,
##     cl = "splice"
##   )
## }
## #' @importFrom grid convertX convertY polylineGrob gList setChildren
## #' @export
## makeContent.splice <- function(x) {
##   data <- x$data
  
##   #data$x <- convertX(data$x, "mm", valueOnly = TRUE)
##   #data$xend <- convertX(data$xend, "mm", valueOnly = TRUE)
##   #data$y <- convertY(data$y, "mm", valueOnly = TRUE)
##   #data$height <- convertY(data$height, "mm", valueOnly = TRUE)
  
##   splices <- lapply(seq_along(data$x), function(i) {
##     cbind(
##       create_splice(
##         x = data$x[i],
##         y = data$y[i],
##         xend = data$xend[i],
##         height = data$height[i],
##         spline_shape = x$spline_shape
##       ),
##       id = i
##     )
##   })
##   splices <- do.call(rbind, splices)
##   splice_paths <- grid::polylineGrob(
##     x = splices$x,
##     y = splices$y,
##     id = splices$id,
##     #default.units = "mm",
##     gp = x$gp
##   )
##   setChildren(x, gList(splice_paths))
## }
