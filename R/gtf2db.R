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

#' @title isGTF
#' @description Check if an object is a GTF database.
#' @param gtf An object to check.
#' @return Logical, TRUE if the object is a GTF database.
#' @export
isGTF <- function(gtf = NULL) {
  return("GTF" %in% class(gtf))
}

#' @title notGTF
#' @description Check if an object is NOT a GTF database.
#' @param gtf An object to check.
#' @return Logical, TRUE if the object is not a GTF database.
#' @export
notGTF <- function(gtf = NULL) {
  return("GTF" %ni% class(gtf))
}
