#' @importFrom ComplexHeatmap Heatmap rowAnnotation HeatmapAnnotation 
#' @export
hclustCells <- function(object = NULL, assay = NULL, layer = "data", features = NULL, cells = NULL, cell.group = 5, feature.group = 5, anno.col.tags = NULL, anno.row.tags = NULL, col.cols = NULL, row.cols = NULL, distance.method = c("euclidean","maximum", "manhattan", "binary", "canberra", "minkowski"), hclust.method = c("complete", "ward.D","ward.D2","single", "average", "mcquitty", "median", "centroid"), binarize = FALSE, cutoff.bin = 0, use_raster = TRUE, scale = c("none", "row", "column"), return.data = FALSE, ...)
{
  assay <- assay %||% DefaultAssay(object)
  data <- GetAssayData1(object, assay = assay, layer = layer)
  features <- features %||% rownames(data)
  cells <- cells %||% colnames(data)

  features <- intersect(rownames(data), features)
  cells <- intersect(colnames(data), cells)

  data <- data[features, cells]

  features <- names(which(rowSums(data>0)!=0))
  cells <- names(which(colSums(data>0)!=0))
  data <- data[features, cells]
   
  if (isTRUE(binarize)) {
    data <- data > cutoff.bin
  }
  data <- as.matrix(data) + 0
 
  scale <- match.arg(scale)
  if (scale == "row") {
    data <- apply(data, 1, base::scale)
    data <- t(data)
    colnames(data) <- cells
  }
  if (scale == "column") {
    data <- apply(data, 2, base::scale)
    rownames(data) <- features
  }

  distance.method <- match.arg(distance.method)
  hclust.method <- match.arg(hclust.method)
  d <- dist(data, method = distance.method)
  hc <- hclust(d, method = hclust.method)
  dend <- as.dendrogram(hc)
  km <- cutree(hc, k = feature.group)
  km <- as.factor(sort(km))
  d <- dist(t(data), method = distance.method)
  hc <- hclust(d, method = hclust.method)
  dend1 <- as.dendrogram(hc)
  km1 <- cutree(hc, k = cell.group)
  km1 <- as.factor(sort(km1))
  
  data <- data[names(km),names(km1)]

  if (isTRUE(return.data)) {
    return(data)
  }
  
  lefta <- rowAnnotation(module = km, border = TRUE, annotation_name_side = "top")
  topa <- HeatmapAnnotation(cluster = km1, border=TRUE)
  righta <- NULL  
  if (!is.null(anno.row.tags)) {
    df <- obj[[assay]][[]]
    tags <- intersect(anno.row.tags, colnames(df))
    if (!is.null(tags)) {
      df <- df[features, tags]
      righta <- rowAnnotation(meta = df, border=TRUE)
    }
  }
  
  if (isTRUE(binarize)) {
    ht <- Heatmap(data, top_annotation = topa, left_annotation = lefta, column_split = km1, row_split = km, border=TRUE, cluster_rows = FALSE, cluster_columns = FALSE, use_raster = use_raster, col = c("white", "black"), right_annotation = righta, ...)
  } else {
    ht <- Heatmap(data, top_annotation = topa, left_annotation = lefta, column_split = km1, row_split = km, border=TRUE, cluster_rows = FALSE, cluster_columns = FALSE, use_raster = use_raster, right_annotation = righta, ...)
  }
  ht
}








