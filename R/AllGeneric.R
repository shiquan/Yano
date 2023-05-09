setGeneric("QuickRecipe", function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                   nvar = 2000, ...) standardGeneric("QuickRecipe"))

setGeneric("QuickRecipe0", function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                    ...) standardGeneric("QuickRecipe0"))

setGeneric("SMerge", function(A = NULL, B = NULL, ..., intersect = c("none", "row", "column", "both")) standardGeneric("SMerge"))

setMethod("SMerge", signature(A = "MatrixOrNull", B = "MatrixOrNull"),function(A, B, ..., intersect = "none") {
  if (is.null(B)) return(A)
  if (is.null(A)) stop("No matrix.")
  
  if (intersect %in% c("row", "both")) {
    row.nm <- intersect(rownames(A), rownames(B))
    A <- A[row.nm,]
    B <- B[row.nm,]
  }
  
  if (intersect %in% c("column", "both")) {
    col.nm <- intersect(colnames(A), colnames(B))
    A <- A[,col.nm]
    B <- B[,col.nm]
  }
  
  row.nm <- unique(c(rownames(A), rownames(B)))
  col.nm <- unique(c(colnames(A), colnames(B)))
  
  A <- as(A, "TsparseMatrix")
  B <- as(B, "TsparseMatrix")

  C <- Matrix::sparseMatrix(i = c(match(rownames(A), row.nm)[A@i+1], match(rownames(B), row.nm)[B@i+1]),
                            j = c(match(colnames(A), col.nm)[A@j+1], match(colnames(B), col.nm)[B@j+1]),
                            x = c(A@x, B@x), dims = c(length(row.nm), length(col.nm)))
  rownames(C) <- row.nm
  colnames(C) <- col.nm

  SMerge(C, ..., intersect=intersect)
} )

setMethod("SMerge", signature(A = "list"),
          definition = function(A = NULL, intersect = "none") {
            A[['intersect']] <- intersect
            do.call(SMerge, A)
          })
