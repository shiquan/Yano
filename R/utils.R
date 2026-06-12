#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @title Meta
#' @description Access the meta feature table of a Seurat assay.
#' @param object Seurat object.
#' @param assay Assay name. If NULL, uses the default assay.
#' @return A data.frame of meta features.
#' @export
Meta <- function(object, assay = NULL)
{
  assay <- assay %||% DefaultAssay(object)
  object[[assay]][[]]
}

# Helper: find default dimensionality reduction
DefaultDimReduc <- function(object) {
  reductions <- Seurat::Reductions(object)
  for (r in c("umap", "tsne", "pca")) {
    if (r %in% reductions) return(r)
  }
  if (length(reductions) > 0) return(reductions[1])
  stop("No dimensionality reduction found.")
}

# Helper: generate random name string
RandomName <- function(length = 10L) {
  paste0("RandomName_", paste(sample(c(letters, 0:9), length, replace = TRUE), collapse = ""))
}

# Helper: Seurat-compatible ggplot theme helpers
# Return rhs when lhs is not NULL, otherwise NULL.
# Semantics: %iff% acts as a guard — only evaluate rhs if lhs is present.
`%iff%` <- function(lhs, rhs) {
  if (!is.null(lhs)) return(rhs)
  return(lhs)
}

CenterTitle <- function() {
  theme(plot.title = element_text(hjust = 0.5))
}

NoLegend <- function() {
  theme(legend.position = "none")
}

NoAxes <- function() {
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )
}

# ---- Tree helpers (replace Seurat:::DFT, MapVals, GetLeftDescendants) ----

# Depth-first traversal: return all descendant node indices from a phylo tree
DFT <- function(tree, node, include.children = TRUE) {
  if (!include.children) return(node)
  children <- tree$edge[tree$edge[, 1] == node, 2]
  if (length(children) == 0) return(node)
  unique(c(node, unlist(lapply(children, function(x) DFT(tree, x, include.children = TRUE)))))
}

# Map values: for each entry in vec, return the corresponding entry in to
# using from as the lookup key (i.e. to[match(vec, from)])
MapVals <- function(vec, from, to) {
  to[match(vec, from)]
}

# Get leftmost descendant tip indices for a given internal node
GetLeftDescendants <- function(tree, node) {
  current <- node
  while (TRUE) {
    children <- tree$edge[tree$edge[, 1] == current, 2]
    if (length(children) == 0) break
    current <- min(children)
  }
  # Return tip index (leaves are 1..Ntip, internal nodes are > Ntip)
  # The leftmost tip is a leaf → return its index in tip.label
  which(tree$tip.label %in% tree$tip.label[current])
}

#' @useDynLib Yano, .registration = TRUE
NULL
