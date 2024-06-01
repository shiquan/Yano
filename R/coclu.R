#' @importFrom Seurat AggregateExpression NormalizeData ScaleData FetchData
#' @export
AggAndCoClust <- function(object = NULL, assay = "EPT", features = NULL,
                         bind.assay = "RNA", bind.name = "gene_name")
{
  agg <- AggregateExpression(object, assays = assay, return.seurat = TRUE)
  agg1 <- AggregateExpression(object, assays = bind.assay, return.seurat = TRUE)

  agg <- NormalizeData(agg)
  agg <- ScaleData(agg, features = features)

  object[[assay]]@meta.features[features,] %>% pull(bind.name) -> bind.features
  
  agg1 <- NormalizeData(agg1)
  agg1 <- ScaleData(agg1, features = unique(bind.features))

  dat <- FetchData(agg, layer = "scale.data", vars = features)
  dat <- t(dat)

  dat1 <- FetchData(agg1, layer = "scale.data", vars = unique(bind.features))
  dat1 <- t(dat1)

  dat <- dat[features,]
  dat1 <- dat1[bind.features,]

  dat0 <- cbind(dat, dat1)
  d <- dist(dat0)
  hc <- hclust(d)
  dat <- dat[hc$order,]
  dat1 <- dat1[hc$order,]

  return(list(dat, dat1, hc))
}

