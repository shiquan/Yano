setGeneric(
  name = "mergeMatrix",
  def=function(x,  ...)  standardGeneric("mergeMatrix")
)

#' @export
setMethod(f = "mergeMatrix",
          signature = signature(x="SMatrix"),
          definition = function(x = NULL, y = NULL, ...) {
            if (is.null(x)) stop("X is empty.")
            if (is.null(y)) return(x)
            
            rnames <- unique(c(rownames(x), rownames(y)))
            cnames <- unique(c(colnames(x), colnames(y)))
            
            rl <- length(rnames)
            cl <- length(cnames)
            
            ra <- c(1:rl)
            ca <- c(1:cl)
            
            names(ra) <- rnames
            names(ca) <- cnames
            
            x <- as(x, "dgTMatrix")
            y <- as(y, "dgTMatrix")
            
            cat(paste0(rl, " features X ", cl, " vars.\n"))
            
            i <- c(ra[rownames(x)[x@i+1]],ra[rownames(y)[y@i+1]])
            j <- c(ca[colnames(x)[x@j+1]],ca[colnames(y)[y@j+1]])
            m <- Matrix::sparseMatrix(i = i, j = j, x = c(x@x, y@x), dims = c(rl, cl))
            rownames(m) <- rnames
            colnames(m) <- cnames
            return(mergeMatrix(m, ...))
          }
          )
