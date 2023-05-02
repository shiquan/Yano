FindAllAltEPT <- function(object, group.by = NULL, subset.ident=NULL,
                          slot = "counts", assay="EPT", block.name = "gene_name",
                          type =  c("exon", "exonintron", "multiexon", "utr3", "utr5"),
                          features=NULL, logfc.threshold = 0.25, test.use = "DEXSeq",
                          min.pct = 0.05, max.cells.per.ident = 100000,
                          min.cells.per.ident = 20,
                          min.cells.per.group = 20,
                          min.cells.feature =3,
                          pesudo.group = 6, random.seed = 1,
                          min.epts.per.block = 2,
                          return.thresh = 1e-2, verbose = TRUE, ...
                          )
{
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }

  idents.all <- sort(x=unique(Idents(object=object)))

  epts.alt <- list()
  messages <- list()
  
  for (i in 1:length(x=idents.all)) {
    if (verbose) {
      message("Caculating alternative EPT expression for cluster ", idents.all[i])
    }

    tab <- FindAltEPT(
      object = object,
      ident.1 = idents.all[[i]],
      ident.2 = NULL,
      group.by = group.by,
      subset.ident = subset.ident,
      features=features,
      assay=assay,
      slot=slot,
      block.name = block.name,
      logfc.threshold = logfc.threshold,
      test.use = test.use,
      min.pct=min.pct,
      type=type,
      max.cells.per.ident=max.cells.per.ident,
      min.cells.per.group = min.cells.per.group,
      min.cells.feature=min.cells.feature,
      pesudo.group = pesudo.group,
      verbose=verbose,
      min.epts.per.block = min.epts.per.block,
      return.thresh = return.thresh,
      ...
    )

    tab$cluster <- i
    epts.alt[[i]] <- tab
    #epts.alt[[i]]
    #epts.alt[[i]]$cluster <- i
  }
  
  return(rbindlist(epts.alt))
}

FindAltEPT.0 <- function(object, cells.1 = NULL, cells.2 = NULL, features = NULL, 
                         slot = "counts", assay = "EPT", block.matrix = NULL, block.name= "gene_name",
                         logfc.threshold = 0.25, test.use = "glm.nb", min.pct = 0.05, verbose = TRUE,
                         max.cells.per.ident = 10000, min.cells.per.ident = 20,
                         min.cells.per.group = 20,
                         type =  c("exon", "exonintron", "multiexon", "utr3", "utr5"),
                         min.cells.feature = 3, pesudo.group = 6, return.thresh = 1e-2,
                         min.epts.per.block = 2,
                         ... )
{

  DefaultAssay(object) <- assay
  
  feature.tab <- object[[assay]]@meta.features

  if (is.null(feature.tab) | length(intersect(c("name", "type"), colnames(feature.tab))) != 2)
    stop("Load annotated bed file first.")

  if (!(block.name %in% colnames(feature.tab)))
    stop("block.name not found in feature table. Check block.name and @meta.feature first.")

  if (is.null(cells.2)) {
    cells.2 <- setdiff(colnames(object), cells.1)
    if (length(cells.2) > max.cells.per.ident) {
      cells.2 <- sample(cells.2, size = max.cells.per.ident)
    }
  }
  Seurat:::ValidateCellGroups(object=object, cells.1=cells.1, cells.2=cells.2,
                              min.cells.group = min.cells.per.ident)

  features <- features %||% rownames(x = object)
  
  n1 <- length(cells.1)
  n2 <- length(cells.2)
  if (n1 < min.cells.per.ident) stop("Too few cells in this ident 1.")
  if (n2 < min.cells.per.ident) stop("Too few cells in this ident 2.")
  
  feature.tab <- feature.tab[,c("name", block.name, "type")]
  feature.tab <- feature.tab[which(feature.tab$type %in% type),]
  tbl <- table(feature.tab[,block.name])
  genes <- names(tbl)[which(tbl >= min.epts.per.block)]
  feature.tab <- feature.tab[which(feature.tab$gene_name %in% genes),]
  features<- intersect(features,feature.tab$name)

  if (is.null(features)) stop("No valid features.")
  ## check feature coverage
  cnt <- Seurat::GetAssayData(object = object, slot = slot, assay= assay)[features,cells.1]
  
  features <- features[which(Matrix::rowSums(cnt>0)/n1 > min.pct)]
  
  feature.tab <- feature.tab[features,]

  group1 <- as.integer(n1 / min.cells.per.ident)
  if (group1 > pesudo.group) group1 <- pesudo.group

  group2 <- as.integer(n2 / min.cells.per.ident)
  if (group2 > pesudo.group) group2 <- pesudo.group
  
  new.group1 <- rep(1:group1,(n1/group1+1), length.out = n1)
  new.group2 <- rep((group1+1):(group1+group2),(n2/group2+1), length.out = n2)

  cnt <- Seurat::GetAssayData(object = object, slot = slot, assay= assay)[features,c(cells.1, cells.2)]

  new.group <- c(new.group1,new.group2)
  cnt <- as(cnt,"TsparseMatrix")
  cnt <- Matrix::sparseMatrix(i=(cnt@i+1),j=new.group[cnt@j+1], x=cnt@x)
  colnames(cnt) <- paste0("group_",1:(group1+group2))
  rownames(cnt)<- features

  meta.tab <- data.frame(row.names=colnames(cnt),
                         condition=c(rep("target", group1),rep("comparison", group2)))

  dxd = DEXSeq::DEXSeqDataSet(as.matrix(cnt),
                              sampleData=meta.tab,
                              groupID = feature.tab$gene_name,
                              featureID = feature.tab$name,
                              design= ~sample+exon+condition:exon)

  dxd = DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
  dxd = DEXSeq::estimateDispersions(dxd)
  dxd = DEXSeq::testForDEU(dxd)
  dxd = DEXSeq::estimateExonFoldChanges(dxd)
  results = DEXSeq::DEXSeqResults(dxd)

  results<-results[order(results$padj),]
  results <- results[which(abs(results$log2fold_target_comparison) > logfc.threshold & results$padj < return.thresh),]

  results <- as.data.frame(results)[, c("groupID", "featureID", "exonBaseMean", "pvalue", "padj", "log2fold_target_comparison")]
  return(results)
}

FindAltEPT <- function(object, ident.1 = NULL, ident.2 = NULL,
                       group.by = NULL,
                       subset.ident = NULL,
                       features = NULL,
                       assay = "EPT", slot = "counts", block.name = "gene_name",
                       logfc.threshold = 0.25, test.use = "DEXSeq", min.pct = 0.1, verbose = TRUE,
                       max.cells.per.ident = 10000, min.cells.per.ident = 20,
                       min.cells.per.group = 20,
                       type =  c("exon", "exonintron", "multiexon", "utr3", "utr5"),
                       min.cells.feature = 3, pesudo.group = 6,
                       return.thresh = 1e-2, min.epts.per.block = 2,
                       ... )
{
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }

  if (length(x = ident.1) == 0) {
    stop("At least 1 ident must be specified in `ident.1`")
  }
  
  cells <- Seurat:::IdentsToCells(object=object, ident.1 = ident.1, ident.2 = ident.2,
                         cellnames.use = colnames(object))
  
  
  alts <- FindAltEPT.0(object=object, cells.1=cells$cells.1, cells.2=cells$cellls.2, features=features,
                       slot=slot, assay=assay, type=type, pesudo.group=pesudo.group,
                       logfc.threshold = logfc.threshold, test.use=test.use, min.pct = min.pct,
                       verbose=verbose, max.cells.per.ident = max.cells.per.ident,
                       min.cells.per.ident = min.cells.per.ident, min.cells.per.group =min.cells.per.group,
                       min.epts.per.block=min.epts.per.block, ...)


  return(alts)
}

