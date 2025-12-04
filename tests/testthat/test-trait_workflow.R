test_that("build_trait_se builds a SummarizedExperiment from toy data", {
  path <- system.file("extdata", "pGlyco3_gpsm_toyexample.txt",
    package = "glycoTraitR"
  )
  gpsm_toyexample <- read_pGlyco3_gpsm(path)

  data("meta_toyexample", package = "glycoTraitR")

  trait_se <- build_trait_se(
    gpsm   = gpsm_toyexample,
    from   = "pGlyco3",
    motifs = NULL,
    level  = "protein",
    meta   = meta_toyexample
  )

  expect_s4_class(trait_se, "SummarizedExperiment")
  expect_gt(nrow(trait_se), 0)
  expect_gt(ncol(trait_se), 0)
  expect_gt(length(SummarizedExperiment::assayNames(trait_se)), 0)
})

test_that("analyze_trait_changes runs on trait_se without errors", {
  path <- system.file("extdata", "pGlyco3_gpsm_toyexample.txt",
    package = "glycoTraitR"
  )
  gpsm_toyexample <- read_pGlyco3_gpsm(path)
  data("meta_toyexample", package = "glycoTraitR")

  trait_se <- build_trait_se(
    gpsm   = gpsm_toyexample,
    from   = "pGlyco3",
    motifs = NULL,
    level  = "protein",
    meta   = meta_toyexample
  )

  changed_traits <- analyze_trait_changes(
    trait_se     = trait_se,
    group_col    = "Diagnosis",
    group_levels = c("Normal", "Symptomatic"),
    min_psm      = 20
  )

  expect_s3_class(changed_traits, "data.frame")
  expect_gt(nrow(changed_traits), 0)
  expect_true(all(c("trait", "level") %in% colnames(changed_traits)))
})

test_that("plot_trait_distribution returns ggplot objects", {
  path <- system.file("extdata", "pGlyco3_gpsm_toyexample.txt",
    package = "glycoTraitR"
  )
  gpsm_toyexample <- read_pGlyco3_gpsm(path)
  data("meta_toyexample", package = "glycoTraitR")

  trait_se <- build_trait_se(
    gpsm   = gpsm_toyexample,
    from   = "pGlyco3",
    motifs = NULL,
    level  = "protein",
    meta   = meta_toyexample
  )

  changed_traits <- analyze_trait_changes(
    trait_se     = trait_se,
    group_col    = "Diagnosis",
    group_levels = c("Normal", "Symptomatic"),
    min_psm      = 20
  )

  trait_name <- changed_traits$trait[1]
  feature <- changed_traits$level[1]

  p <- plot_trait_distribution(
    trait_se     = trait_se,
    group_col    = "Diagnosis",
    group_levels = c("Normal", "Symptomatic"),
    trait_name   = trait_name,
    feature      = feature
  )

  # plot_trait_distribution returns a list with two ggplots: p_hist and p_box
  expect_true(is.list(p))
  expect_true(all(c("p_hist", "p_box") %in% names(p)))
  expect_s3_class(p$p_hist, "ggplot")
  expect_s3_class(p$p_box, "ggplot")
})
