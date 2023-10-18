#' @export matrix.dist
matrix.dist <- function(matrix = NULL, method = "Jaccard")
{
  d <- .Call("matrix_distance", matrix, method);
  return(d)
}
