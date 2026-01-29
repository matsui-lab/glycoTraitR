test_that("build_trait_se produces a SummarizedExperiment", {
  path <- system.file("extdata", "pGlyco3_gpsm_toyexample.txt",
    package = "glycoTraitR"
  )

  gpsm <- read_pGlyco3_gpsm(path)
  data("meta_toyexample", package = "glycoTraitR")

  se <- build_trait_se(
    gpsm   = gpsm,
    from   = "pGlyco3",
    motifs = NULL,
    level  = "protein",
    meta   = meta_toyexample
  )

  expect_s4_class(se, "SummarizedExperiment")
  expect_true(length(SummarizedExperiment::assays(se)) > 0) # at least one trait
  expect_true(nrow(se) > 0) # at least one protein
})
