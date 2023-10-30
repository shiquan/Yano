# Yano
Yano represents an R/C toolkit designed for conducting spatial dissimilarity analysis on single-cell RNA sequencing data. This method revolves around the core concept of examining the distinct expression patterns of a given feature (e.g. exon, EPT or EAT) in relation to its associated binding feature (typically a gene or genomic locus) within the context of cell lineage (1D), spatial position (2D), or the multi-dimensional PCA space. The discernible differences in feature expression patterns and their binding features provide insights into a range of biological phenomena, including alternative splicing, cis-antisense RNA regulation, allele-specific gene expression, and more.

Given the inherent sparsity and heterogeneity of single-cell RNA data, a precise and informative measure of dissimilarity between two features becomes essential. In many cases, the spatial autocorrelation is a notable feature of non-random expression patterns across cell lineage, spatial location, and UMAP projections. To address this, Yano introduces a novel statistical metric known as the "D score" which combines the spatial autocorrelation of feature X with the dissimilarity between feature X and its associated binding feature Y. The resulting D score follows a normal distribution, allowing for the calculation of a p-value for each X-Y pair with a permutation method (see manuscript for details). 

Yano is seamlessly integrated with Seurat, building upon the Seurat object's framework. Users can perform conventional cell clustering analyses using the state-of-the-art Seurat pipeline and then incorporate exon, EPT, or EAT counts as new "assays" within the Seurat objects. Subsequently, Yano facilitates the assessment of spatial dissimilarity between these two assays. The whole pipeline list below.

![pipeline](https://github.com/shiquan/Yano-doc/blob/master/figs/pipeline.png?raw=true)

# INSTALL

```
devtools::install_package("shiquan/Yano")
```

# Quick start

# Real cases
* Case 1, EXON/EPT/EAT analysis for Human Brain 3k cells
* Case 2, EXON/EPT analysis for 10X Visium mouse brain
* Case 3, EXON/EPT/EAT analysis for multiple samples
* Case 4, EXON analysis for human brain (Smartseq2)

# Issues report

# FAQ

# LICENSE
MIT.

# Citation
 


