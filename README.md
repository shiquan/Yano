# Yano
Yano is an R/C package to perform spatial dissimilarity analysis for single cell RNA sequencing data. The key idea of this method is to test the different expression pattern of a feature (such as exon/EPT/EAT) with its binding feature (gene or locus) in cell lineage (1D), spatial location (2D), or PCA space (high D). The different expressed pattern of a feature and its binding feature indicate various biological phenomenon, such as alternative splcing, cis-antisense RNA regulation, allele-specific gene expression and so on. 

![pipeline](https://github.com/shiquan/Yano-doc/blob/master/figs/pipeline.png?raw=true)

# INSTALL

```
devtools::install_package("shiquan/Yano")
```
To call EPTs and EATs from BAM file, you should also install [PISA](www.github.com/shiquan/PISA).
```
git clone www.github.com/shiquan/PISA
cd PISA
make
```
# Demo
* Case 1, EXON/EPT/EAT analysis for Human Brain 3k cells
* Case 2, EXON/EPT analysis for 10X Visium mouse brain
* Case 3, EXON/EPT/EAT analysis for human PBMC 10k cells
* Case 4, EXON analysis for human brain (Smartseq2)

# Issues report

# FAQ

# LICENSE
MIT.

# Citation
 


