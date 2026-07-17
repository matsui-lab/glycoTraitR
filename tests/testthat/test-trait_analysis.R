testthat::test_that(
  "analyze_hscore_changes returns consistent output structure",
  {
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
      min_samples = 3
    )

    res_ab_2 <- analyze_hscore_changes(
      gpsm = gpsm,
      from = "pGlyco3",
      motifs = NULL,
      meta = meta,
      group_col = "Diagnosis",
      group_levels = c("Normal", "Symptomatic"),
      B = 10,
      min_samples = 3
    )

    testthat::expect_s3_class(res_ab_1, "data.frame")
    testthat::expect_s3_class(res_ab_2, "data.frame")

    testthat::expect_gt(nrow(res_ab_1), 0L)
    testthat::expect_gt(ncol(res_ab_1), 0L)

    testthat::expect_equal(
      nrow(res_ab_1),
      nrow(res_ab_2)
    )

    testthat::expect_equal(
      ncol(res_ab_1),
      ncol(res_ab_2)
    )

    testthat::expect_identical(
      names(res_ab_1),
      names(res_ab_2)
    )

    column_classes_1 <- vapply(
      res_ab_1,
      function(x) paste(class(x), collapse = "/"),
      character(1)
    )

    column_classes_2 <- vapply(
      res_ab_2,
      function(x) paste(class(x), collapse = "/"),
      character(1)
    )

    testthat::expect_identical(
      column_classes_1,
      column_classes_2
    )

    required_columns <- c(
      "trait",
      "score_type",
      "feature",
      "level",
      "diff",
      "pval"
    )

    testthat::expect_true(
      all(required_columns %in% names(res_ab_1))
    )
  }
)
