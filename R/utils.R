#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @export
Meta <- function(object, assay = NULL)
{
  assay <- assay %||% DefaultAssay(object)
  object[[assay]][[]]
}
  
#' @useDynLib Yano, .registration = TRUE
NULL
