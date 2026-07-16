testthat::test_that("analyze_hscore_changes is reproducible and respects group order", {
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

  old_pb <- pbapply::pboptions(type = "none")
  on.exit(pbapply::pboptions(old_pb), add = TRUE)

  res_ab_1 <- analyze_hscore_changes(
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

  res_ab_2 <- analyze_hscore_changes(
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

  testthat::expect_equal(res_ab_1, res_ab_2)

  res_ba <- analyze_hscore_changes(
    gpsm = gpsm,
    from = "pGlyco3",
    motifs = NULL,
    meta = meta,
    group_col = "Diagnosis",
    group_levels = c("Symptomatic", "Normal"),
    B = 10,
    min_samples = 3,
    seed = 123
  )

  key <- c("trait", "score_type", "feature", "level")

  ab <- res_ab_1[, c(key, "diff")]
  ba <- res_ba[, c(key, "diff")]

  names(ab)[names(ab) == "diff"] <- "diff_ab"
  names(ba)[names(ba) == "diff"] <- "diff_ba"

  paired <- merge(ab, ba, by = key)

  testthat::expect_gt(nrow(paired), 0L)

  testthat::expect_equal(
    paired$diff_ab,
    -paired$diff_ba,
    tolerance = 1e-10
  )
})
