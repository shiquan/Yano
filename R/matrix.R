setGeneric(
  name = "mergeMatrix",
  def=function(x,  ...)  standardGeneric("mergeMatrix")
)

MMerge0 <- function(x = NULL, y = NULL, ...) {
  if (is.null(x)) stop("X is empty.")
  if (is.null(y)) return(x)
  if (is.list(x)) {
    l <- length(x)
    x[[l+1]] <- as(y, "dgCMatrix")
  } else {
    x <- list(as(x, "dgCMatrix"), as(y, "dgCMatrix"))
  }
  MMerge0(x, ...)
}
 
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
