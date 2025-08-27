#' @export
GetAssayData1 <- function(object, assay = NULL, layer = "counts")
{
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay    
  if (packageVersion("Seurat") < numeric_version(as.character(5))) {
    data <- GetAssayData(object, slot = layer)
  } else {
    data <- NULL
    layers <- Layers(object = object[[assay]], search = layer)
    for (i in seq_along(along.with = layers)) {
      l <- layers[i]
      data0 <- LayerData(object[[assay]], layer = l)
      if ("RenameDims" %in% class(data0)) {
        data0 <- as.sparse(data0)
      }
      data <- mergeMatrix(data0, data)
    }
  }
  DefaultAssay(object) <- old.assay
  data
}

FetchData1 <- function(object, assay = NULL, layer, vars = NULL, ...)
{
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay
  if (packageVersion("Seurat") < numeric_version(as.character(5))) {
    dt <- Seurat::FetchData(object, slot = layer, vars = vars, ...)
  } else {
    dt <- Seurat::FetchData(object, layer = layer, vars = vars, ...)
  }

  DefaultAssay(object) <- old.assay
  dt
}
