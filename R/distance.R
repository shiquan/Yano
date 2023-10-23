#' @export matrix.dist
matrix.dist <- function(matrix = NULL, method = "Jaccard", threads=1)
{
  d <- .Call("matrix_distance", matrix, method, threads);  
  return(d)
}
