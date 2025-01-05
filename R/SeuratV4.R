GetAssayData1 <- function(object, assay = NULL, layer = "counts", ...)
{
  assay <- assay %||% DefaultAssay(object)
  old.assay <- DefaultAssay(object)
  DefaultAssay(object) <- assay    
  if (packageVersion("Seurat") < numeric_version(as.character(5))) {
    data <- GetAssayData(object, slot = layer, ...)
  } else {
    data <- NULL
    layers <- Layers(object = object[[assay]], search = layer)

    for (i in seq_along(along.with = layers)) {
      l <- layers[i]
      data0 <- LayerData(object, layer = l)
      data <- mergeMatrix(data0, data)
    }
  }
  DefaultAssay(object) <- old.assay
  data
}

FetchData1 <- function(object, layer, ...)
{
  if (packageVersion("Seurat") < numeric_version(as.character(5))) {
    Seurat::FetchData(object, slot = layer, ...)
  } else {
    Seurat::FetchData(object, layer = layer, ...)
  }  
}
