pkgname <- "Yano"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Yano')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("FindAllAltExp")
### * FindAllAltExp

flush(stderr()); flush(stdout())

### Name: FindAllAltExp
### Title: Test alternative expression for all cell groups
### Aliases: FindAllAltExp

### ** Examples

data("glbt_small")
DefaultAssay(glbt_small) <- "exon"
alt.exon <- FindAllAltExp(object = glbt_small, bind.assay = "RNA", bind.name = "gene_name", features = rownames(glbt_small))
head(alt.exon)




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
