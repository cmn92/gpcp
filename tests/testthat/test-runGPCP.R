test_that("runGPCP works as expected (diploid VCF/HapMap)", {
  phenotypeFile <- read.csv("~/gpcp/data/phenotypeFile.csv")
  genotypeFile  <- system.file("extdata", "genotypeFile_Chr9and11.vcf",
                               package = "gpcp")

  result <- suppressWarnings(
    runGPCP(
      phenotypeFile = phenotypeFile,
      genotypeFile   = genotypeFile,
      genotypes      = "Accession",
      traits         = c("YIELD", "DMC"),
      weights        = c(3, 1),
      userFixed      = c("LOC", "REP"),
      Ploidy         = 2,
      NCrosses       = 150
    )
  )

  expect_s3_class(result, "data.frame")
  expect_true(ncol(result) >= 3)        # at least Parent1, Parent2, Mean, etc.
  expect_equal(nrow(result), 150)       # we asked for top 150 crosses
})

test_that("runGPCP handles a 4-ploid test dataset", {
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

  # 5) expectations
  expect_s3_class(result4x, "data.frame")
  expect_equal(nrow(result4x), 100)
  expect_true(ncol(result4x) >= 3)        # at least Parent1, Parent2, Mean, etc.
})
