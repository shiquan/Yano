GetAssayData1 <- function(object, assay = NULL, layer = "counts", ...)
{
  assay <- assay %||% DefaultAssay(object)
  if (packageVersion("Seurat") < numeric_version(as.character(5))) {
    Seurat::GetAssayData(object, assay = assay, slot = layer, ...)
  } else {
    Seurat::GetAssayData(object, assay = assay, layer = layer, ...)
  }
}

FetchData1 <- function(object, layer, ...)
{
  if (packageVersion("Seurat") < numeric_version(as.character(5))) {
    Seurat::FetchData(object, slot = layer, ...)
  } else {
    Seurat::FetchData(object,layer = layer, ...)
  }  
}
