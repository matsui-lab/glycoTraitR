testthat::test_that("analyze_trait_changes", {
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  testthat::skip_if_not_installed("pbapply")
  testthat::skip_if_not_installed("car")

  # silence pbapply progress bar in tests
  if (requireNamespace("pbapply", quietly = TRUE)) {
    old_pb <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(old_pb), add = TRUE)
  }

  # -----------------------------
  # Construct a minimal SE object
  # -----------------------------
  # 6 samples, 2 groups (A/B), 2 levels (rows), 2 traits (assays)
  samples <- paste0("s", 1:6)
  grp <- factor(c("A", "A", "A", "B", "B", "B"), levels = c("A", "B"))
  coldata <- S4Vectors::DataFrame(grp = grp, row.names = samples)

  rowdata <- S4Vectors::DataFrame(level = c("L1", "L2"), row.names = c("L1", "L2"))

  # Trait 1 matrix
  # Row L1: strong group difference, includes one NA to hit NA-filter path
  # Row L2: all zeros -> should be skipped by all(value_x == 0)
  trait1 <- matrix(
    c(
      10, 11, NA,  1,  2,  3,   # L1
      0,  0,  0,  0,  0,  0    # L2 (all-zero)
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("L1", "L2"), samples)
  )

  # Trait 2 matrix
  # Row L1: similar distributions -> not significant -> should return NULL (else branch)
  # Row L2: group B has < min_psm non-NA -> should be skipped by min_psm check
  trait2 <- matrix(
    c(
      5, 5, 6, 5, 6, 5,          # L1 (non-significant but non-constant)
      1, 1, 1, NA, NA, 2         # L2 (B has only 1 non-NA if min_psm=2)
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("L1", "L2"), samples)
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = S4Vectors::SimpleList(trait1 = trait1, trait2 = trait2),
    rowData = rowdata,
    colData = coldata
  )

  # ---------------------------------------
  # 1) group_col missing -> stop() branch
  # ---------------------------------------
  testthat::expect_error(
    analyze_trait_changes(se, group_col = "BADCOL", group_levels = c("A", "B"), min_psm = 2),
    "was not found in colData"
  )

  # ----------------------------------------------------
  # 2) main run: hits all internal return(NULL) branches
  # ----------------------------------------------------
  res <- analyze_trait_changes(se, group_col = "grp", group_levels = c("A", "B"), min_psm = 2)

  # Expect at least one significant result (trait1 at L1)
  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_true(nrow(res) >= 1)

  # check required columns
  testthat::expect_true(all(c("trait", "level", "l_pval", "f_val", "t_pval", "t_val") %in% names(res)))

  # The only expected significant row is trait1 at L1 (others are skipped or non-significant)
  testthat::expect_true(any(res$trait == "trait1" & res$level == "L1"))

  # sanity: p-values are numeric and finite for returned rows
  testthat::expect_type(res$l_pval, "double")
  testthat::expect_type(res$t_pval, "double")
  testthat::expect_true(all(is.finite(res$l_pval)))
  testthat::expect_true(all(is.finite(res$t_pval)))
})
