#' @title ImputationByWeight
#' @param X Expression counts.
#' @param cells Cells to imputate.
#' @param W Weight matrix.
#' @export
ImputationByWeight <- function(X = NULL, cells = NULL, W = NULL, filter = 0.1)
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

  wcells <- colnames(W)
  cells0 <- colnames(X)
  cells1 <- intersect(cells0, wcells)
  if (length(cells1) != length(wcells)) {
    stop("Weight matrix has cells not in X.")
  }
  if (length(cells0) > length(cells1)) {
    X <- X[,cells1]
  }
  
  if (is.null(cells)) {
    cells <- colnames(X)
  }

  idx <- match(cells, colnames(X))
  names(idx) <- cells
  idx <- sort(idx)
  new.cells <- names(idx)
  X0 <- .Call("imputation1", X, idx, W, filter)
  rownames(X0) <- rownames(X)
  colnames(X0) <- new.cells
  X0[,cells]
}
