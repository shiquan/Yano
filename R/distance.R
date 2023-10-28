#' @export MatrixDist
MatrixDist <- function(matrix = NULL, method = "Jaccard", threads=4)
{
  d <- .Call("matrix_distance", matrix, method, threads);
  colnames(d) <- colnames(matrix)
  rownames(d) <- colnames(matrix)
  return(d)
}
