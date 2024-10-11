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

- [Yano's user guide](https://shiquan.github.io/Yano.html)
- [Report issues](https://github.com/shiquan/Yano/issues)
- [Discussions](https://github.com/shiquan/Yano/discussions)

# Short cases

- [Alternative splicing analysis for scRNA-seq](https://shiquan.github.io/Yano_AS.html)
- [Allele-specific gene expression analysis for scRNA-seq](https://shiquan.github.io/Yano_ASE.html)
- [Annotating and prioritizing genetic variants for scRNA-seq](https://shiquan.github.io/Yano_anno.html)


# Citation
 


