setGeneric(
  name = "mergeMatrix",
  def=function(x,  ...)  standardGeneric("mergeMatrix")
)

MMerge0 <- function(x = NULL, y = NULL, ...) {
  if (is.null(x)) stop("X is empty.")
  if (is.null(y)) {
    return(x)
  }
  if (is.list(x)) {
    l <- length(x)
    x[[l+1]] <- as(y, "dgCMatrix")
  } else {
    x <- list(as(x, "dgCMatrix"), as(y, "dgCMatrix"))
  }
  MMerge0(x, ...)
}
#' @title mergeMatrix
#' @description Merge multiple matrix files into one. At least two matrix files should be specified. Records with the same row name and column name will be summed up.
#' @param x Matrix 1 or a list of Matrix
#' @param y Matrix 2
#' @param ... More matrix files.
#' @export
setMethod(f = "mergeMatrix",
          signature = signature(x="SMatrix"),
          definition = function(x = NULL, y = NULL, ...) {
            if (is.null(x)) stop("X is empty.")
            if (is.null(y)) return(x)

            mlst <- MMerge0(x, y, ...)
            O <- .Call("merge_matrix", mlst)
            as(O, "CsparseMatrix")
          }
          )

setMethod(f = "mergeMatrix",
          signature = signature(x="list"),
          definition = function(x = NULL, ...) {
            if (is.null(x)) stop("X is empty.")
            O <- .Call("merge_matrix", mlst)
            as(O, "CsparseMatrix")
          }
          )
