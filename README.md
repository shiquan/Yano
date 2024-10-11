# Yano
Yano represents an R/C toolkit designed for conducting spatial dissimilarity analysis on single-cell RNA sequencing data. This method revolves around the core concept of examining the distinct expression patterns of a given feature (e.g. exon, junction or genetic variants) in relation to its associated binding feature (typically a gene or genomic locus) within the context of cell lineage (1D), spatial position (2D), or the multi-dimensional PCA space. The discernible differences in feature expression patterns and their binding features provide insights into a range of biological phenomena, including alternative splicing, cis-antisense RNA regulation, allele-specific gene expression, and more.

# INSTALL

```
# Dependencies
install.packages(c("Seurat", "R.utils", "viridis", "devtools"), repos = "https://cloud.r-project.org")
```

```
# Install Yano 
devtools::install_github("shiquan/Yano")
```

# Quick start
```{r}
require(Yano)
data("glbt_small")
DefaultAssay(glbt_small) <- "RNA"
glbt_small <- NormalizeData(glbt_small) %>% FindNeighbors() %>% RunUMAP(dim = 1:20)

DimPlot(glbt_small, label = TRUE, label.size = 5)

DefaultAssay(glbt_small) <- "exon"
glbt_small <- NormalizeData(glbt_small)
glbt_small <- ParseExonName(glbt_small)
glbt_small <- RunAutoCorr(glbt_small)
glbt_small <- SetAutoCorrFeatures(glbt_small)
glbt_small <- RunBlockCorr(glbt_small, bind.name = "gene_name", bind.assay = "RNA")

FbtPlot(glbt_small, val = "gene_name.padj")

glbt_small[['exon']][[]] %>% filter(gene_name.padj <1e-10)

FeaturePlot(glbt_small, features = c("chr19:16095264-16095454/+/TPM4", "TPM4"), order=TRUE)

db <- gtf2db("./gencode.v44.annotation.gtf.gz")
TrackPlot(bamfile="./Parent_SC3v3_Human_Glioblastoma_possorted_genome_bam.bam", gtf =db, gene = "TPM4", junc = TRUE, cell.group = Idents(glbt_small), highlights = c(16095264,16095454))
```
# Short cases

- [Alternative splicing analysis for scRNA-seq](https://shiquan.github.io/Yano_AS.html)
- [Allele-specific gene expression analysis for scRNA-seq](https://shiquan.github.io/Yano_ASE.html)
- [Annotating and prioritizing genetic variants for scRNA-seq](https://shiquan.github.io/Yano_anno.html)

# Issues report

# FAQ

# LICENSE
MIT.

# Citation
 


