#'@export
gtf2db <- function(filename = NULL, use_utr = FALSE) {
  if (is.null(filename)) stop("No gtf file.")
  db <- .Call("gtf2db", normalizePath(filename), use_utr)
  class(db) <- "GTF"
  return(db)
}

#'@export
isGTF <- function(db = NULL) {
  return("GTF" %in% class(db))
}

#'@export
notGTF <- function(db = NULL) {
  return("GTF" %ni% class(db))
}
