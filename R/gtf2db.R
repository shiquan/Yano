#'@export
gtf2db <- function(filename = NULL) {
  if (is.null(filename)) stop("No gtf file.")
  db <- .Call("gtf2db", normalizePath(filename))
  class(db) <- "GTF"
  return(db)
}

#'@export
isGTF <- function(db = NULL) {
  return("GTF" %in% class(db))
}
