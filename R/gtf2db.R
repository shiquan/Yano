#' @title gtf2db
#' @description Generate GTF object from gtf file.
#' @param filename Path to gtf file.
#' @param use_utr Load CDS records or not. Default is FALSE.
#' @return A point to GTF struct.
#' @export
gtf2db <- function(filename = NULL, use_utr = FALSE) {
  if (is.null(filename)) stop("No gtf file.")
  db <- .Call("gtf2db", normalizePath(filename), use_utr)
  class(db) <- "GTF"
  return(db)
}

#'@export
isGTF <- function(gtf = NULL) {
  return("GTF" %in% class(gtf))
}

#'@export
notGTF <- function(gtf = NULL) {
  return("GTF" %ni% class(gtf))
}
