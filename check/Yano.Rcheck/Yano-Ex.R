pkgname <- "Yano"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Yano')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("ParseExonName")
### * ParseExonName

flush(stderr()); flush(stdout())

### Name: ParseExonName
### Title: ParseExonName
### Aliases: ParseExonName

### ** Examples

data("glbt_small")
DefaultAssay(glbt_small) <- "exon"
# Check the meta table before parsing
head(glbt_small[['exon']][[]])

glbt_small <- ParseExonName(glbt_small)

# Now see the meta table after parsing
head(glbt_small[['exon']][[]])




cleanEx()
nameEx("RunSDT")
### * RunSDT

flush(stderr()); flush(stdout())

### Name: RunSDT
### Title: RunSDT
### Aliases: RunSDT

### ** Examples

data("glbt_small")
DefaultAssay(glbt_small) <- "RNA"
glbt_small <- NormalizeData(glbt_small) %>% RunUMAP(dim = 1:20)
DefaultAssay(glbt_small) <- "exon"
glbt_small <- NormalizeData(glbt_small)
glbt_small <- ParseExonName(glbt_small)
glbt_small <- RunAutoCorr(glbt_small)
glbt_small <- RunSDT(glbt_small, bind.name = "gene_name", bind.assay = "RNA")




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
