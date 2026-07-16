test_that("analyze_hscore_changes runs successfully", {

  gpsm <- readRDS(system.file(
    "extdata",
    "gpsm_toyexample.rds",
    package = "glycoTraitR"
  ))

  meta <- readRDS(system.file(
    "extdata",
    "meta_toyexample.rds",
    package = "glycoTraitR"
  ))

  res <- analyze_hscore_changes(
    gpsm = gpsm,
    from = "pGlyco3",
    motifs = NULL,
    meta = meta,
    group_col = "Diagnosis",
    group_levels = c("Normal", "Symptomatic"),
    B = 10,
    min_samples = 3,
    seed = 123
  )

  ## returns
  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)

  ## documented output columns
  expect_true(all(c(
    "trait",
    "score_type",
    "feature",
    "diff",
    "pval",
    "level"
  ) %in% names(res)))

  ## p values are valid
  expect_true(all(res$pval >= 0 & res$pval <= 1))

  ## supported levels
  expect_true(all(res$level %in% c("site", "protein")))
})
