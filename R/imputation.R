#' @title ImputationByWeight
#' @param X Expression counts.
#' @param cells Cells to imputate.
#' @param W Weight matrix.
#' @param filter Value below this cutoff will be filtered. Used to reduce density of matrix.
#' @export
ImputationByWeight <- function(X = NULL, cells = NULL, W = NULL, filter = 0.001)
{
  if (is.null(X) || is.null(W)) {
    stop("X and/or W is not set.")
  }

  if ("dgCMatrix" %ni% class(X)) {
    X <- as(X, "CsparseMatrix")
  }
  if ("dgCMatrix" %ni% class(W)) {
    W <- as(x, "CsparseMatrix")
  }
  
  wcells <- rownames(W)
  cells0 <- colnames(X)
  cells1 <- intersect(cells0, wcells)
  if (length(cells1) != length(wcells)) {
    stop("Weight matrix has cells not in X.")
  }

  X <- X[,wcells]
  
  cells <- cells %||% colnames(W)
  cells0 <- intersect(cells, colnames(X))
  idx <- match(cells, colnames(W))
  names(idx) <- cells
  idx <- sort(idx)
  new.cells <- names(idx)

  X0 <- .Call("imputation1", X, idx, W, filter)
  rownames(X0) <- rownames(X)
  colnames(X0) <- new.cells
  
  X0[,cells]
}
