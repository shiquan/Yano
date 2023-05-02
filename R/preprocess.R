#' Quick clust single cell gene expression matrix by Seurat pipeline
#'
#' @rdname QuickRecipe
#' @import Seurat
#' @import Matrix
#'
#' @export
#' 
setMethod(f = "QuickRecipe",
          signature = signature(counts = "SMatrix"),
          definition = function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                nvar = 3000, resolution = 0.5, default.assay = "RNA",
                                new.assay.name = "SCT", ndim=20, ...
                                ) {
            
            seu <- CreateSeuratObject(counts = counts, min.cells = min.cells, min.features = min.features, assay = default.assay)
            QuickRecipe(seu, meta.data = meta.data, min.cells =min.cells, min.features = min.features,
                        nvar = nvar, resolution = resolution, default.assay = default.assay,
                        new.assay.name = new.assay.name, ndim=ndim,...)
          })


setMethod(f = "QuickRecipe",
          signature = signature(counts = "Seurat"),
          definition = function(counts = NULL, meta.data = NULL, min.cells = 20, min.features = 200,
                                nvar = 3000, resolution = 0.5, default.assay = "RNA",
                                new.assay.name = "SCT", ndim = 20,
                                ...
                                ) {
            
            #seu <- CreateSeuratObject(counts = counts, min.cells = min.cells, min.features = min.features)
            seu <- SCTransform(counts, variable.features.n = nvar, assay = default.assay,
                               new.assay.name = new.assay.name)
            seu <- RunPCA(seu, featureres = VariableFeatures(object = seu))
            seu <- FindNeighbors(seu, dims = 1:ndim)
            seu <- FindClusters(seu, resolution = resolution)
            seu <- RunUMAP(seu, dims = 1:ndim)
            seu
          })

