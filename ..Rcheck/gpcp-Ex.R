pkgname <- "gpcp"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "gpcp-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('gpcp')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("phenotypeFile")
### * phenotypeFile

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: phenotypeFile
### Title: Example Phenotype Data
### Aliases: phenotypeFile
### Keywords: datasets

### ** Examples

data(phenotypeFile)
head(phenotypeFile)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("phenotypeFile", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("runGPCP")
### * runGPCP

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: runGPCP
### Title: Genomic Prediction of Cross Performance This function performs
###   genomic prediction of cross performance using genotype and phenotype
###   data.
### Aliases: runGPCP

### ** Examples

# Load phenotype data from CSV
# Diploid pipeline
phenotypeFile <- read.csv(system.file("extdata", "phenotypeFile.csv", package = "gpcp"))
genotypeFile <- system.file("extdata", "genotypeFile_Chr9and11.vcf", package = "gpcp")
finalcrosses <- runGPCP(
    phenotypeFile = phenotypeFile,
    genotypeFile = genotypeFile,
    genotypes = "Accession",
    traits = "YIELD,DMC",
    weights = c(3, 1),
    userFixed = "LOC,REP",
    Ploidy = 2,
    NCrosses = 150
)
message(finalcrosses)
 #PolyPLoid Pipeline
 # 1) load example data from the package
data(DT_polyploid, package = "sommer")
DT <- DT_polyploid
GT <- GT_polyploid
MP <- MP_polyploid
# 2) convert A/T/C/G strings to numeric codes
numo <- sommer::atcg1234(data = GT, ploidy = 4)

# 3) find the set of individuals common to genotypes and phenotypes
common <- intersect(DT$Name, rownames(numo$M))
marks  <- numo$M[common, , drop = FALSE]
pheno2 <- as.data.frame(DT[match(common, DT$Name), ])
# 4) call runGPCP with ploidy = 4
result4x <- suppressWarnings(
  runGPCP(
    phenotypeFile = pheno2,
    genotypeData   = marks,
    genotypes      = "Name",
    traits         = c("total_yield", "tuber_length"),
    weights        = c(3, 1),
    Ploidy         = 4,
    NCrosses       = 100
  )
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("runGPCP", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
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
