\name{runGPCP}
\alias{runGPCP}
\title{
Genomic Prediction of Cross Performance
}
\description{
This function performs genomic prediction of cross performance using genotype and phenotype data. It processes data in several steps including loading necessary software, converting genotype data, processing phenotype data, fitting mixed models, and predicting cross performance based on weighted marker effects.
}
\usage{
runGPCP(phenotypeFile, genotypeFile, genotypes, traits, weights = NA, userSexes = "",
        userFixed = NA, userRandom = NA, Ploidy = NA, NCrosses = NA)
}
\arguments{
  \item{phenotypeFile}{
    A data frame containing phenotypic data, typically read from a CSV file.
  }
  \item{genotypeFile}{
    A file path to the genotypic data, either in VCF format or as a HapMap.
  }
  \item{genotypes}{
    A character string representing the column name in the phenotype file that corresponds to the genotype IDs.
  }
  \item{traits}{
    A string of comma-separated trait names from the phenotype file, which will be used for genomic prediction.
  }
  \item{weights}{
    A numeric vector specifying the weights for the traits. The order of weights should correspond to the order of traits.
  }
  \item{userSexes}{
    Optional. A string representing the column name in the phenotype file corresponding to the individuals' sexes.
  }
  \item{userFixed}{
    A string of comma-separated fixed effect variables from the phenotype file. If no fixed effects are required, set to NA.
  }
  \item{userRandom}{
    A string of comma-separated random effect variables from the phenotype file. If no random effects are required, set to NA.
  }
  \item{Ploidy}{
    An integer representing the ploidy level of the organism (e.g., 2, 4, 6).
  }
  \item{NCrosses}{
    An integer specifying the number of top crosses to output. Maximum is a full diallel.
  }
}
\details{
This function is designed for genomic prediction of cross performance and can handle both diploid and polyploid species. It processes genotype data, calculates genetic relationships, and fits mixed models using the `sommer` package. It outputs the best predicted crosses based on user-defined traits and weights.
}
\value{
A data frame containing predicted crosses with the following columns:
  \item{Parent1}{First parent genotype ID.}
  \item{Parent2}{Second parent genotype ID.}
  \item{CrossPredictedMerit}{Predicted merit of the cross.}
  \item{P1Sex}{Optional. Sex of the first parent if userSexes is provided.}
  \item{P2Sex}{Optional. Sex of the second parent if userSexes is provided.}
}
\references{
Xiang, J., et al. (2016). "Mixed Model Methods for Genomic Prediction." \emph{Nature Genetics}.
Batista, L., et al. (2021). "Genetic Prediction and Relationship Matrices." \emph{Theoretical and Applied Genetics}.
}
\author{
Marlee Labroo, Christine Nyaga, Lukas Mueller
}
\note{
This function relies on the `sommer`, `dplyr`, and `AGHmatrix` packages for processing mixed models and genomic data.
}
\seealso{
\code{\link{sommer}},\code{\link{dplyr}},\link[AGHmatrix]{Gmatrix}
}
\examples{
# Load phenotype data from CSV
phenotypeFile <- read.csv("~/Documents/GCPC_input_files/2020_TDr_PHENO (1).csv")

# Genotype file path
genotypeFile <- "~/Documents/GCPC_input_files/genotypeFile.vcf"


# Define inputs
genotypes <- "Accession"
traits <- c("rAUDPC_YMV", "YIELD", "DMC")
weights <- c(0.2, 3, 1)
userFixed <- c("LOC", "REP")
Ploidy <- 2
NCrosses <- 150

# Run genomic prediction of cross performance
finalcrosses <- runGPCP(
    phenotypeFile = phenotypeFile,
    genotypeFile = genotypeFile,
    genotypes = genotypes,
    traits = paste(traits, collapse = ","),
    weights = weights,
    userFixed = paste(userFixed, collapse = ","),
    Ploidy = Ploidy,
    NCrosses = NCrosses
)

# View the predicted crosses
print(finalcrosses)
}
\keyword{genomic prediction}
\keyword{cross performance}
