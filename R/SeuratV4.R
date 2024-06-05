GetAssayData1 <- function(object, assay = "RNA", layer = "counts", ...)
{
  if (packageVersion("Seurat") < 5) {
    Seurat::GetAssayData(object, assay = assay, slot = layer, ...)
  } else {
    Seurat::GetAssayData(object, assay = assay, layer = layer, ...)
  }
}

FetchData1 <- function(object, layer, ...)
{
  if (packageVersion("Seurat") < 5) {
    Seurat::FetchData(object, slot = layer, ...)
  } else {
    Seurat::FetchData(object,layer = layer, ...)
  }  
}
