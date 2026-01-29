testthat::test_that("plot_trait_distribution returns valid ggplot objects", {
  testthat::skip_if_not_installed("SummarizedExperiment")
  testthat::skip_if_not_installed("ggplot2")

  # -----------------------------
  # Construct a minimal SE object
  # -----------------------------
  # 6 PSMs, 2 groups
  assays <- list(
    TraitA = matrix(
      c(10, 12, NA, 3, 4, 5),
      nrow = 1,
      dimnames = list("Protein1", paste0("psm", 1:6))
    )
  )

  coldata <- S4Vectors::DataFrame(
    Group = c("A", "A", "A", "B", "B", "B"),
    row.names = paste0("psm", 1:6)
  )

  rowdata <- S4Vectors::DataFrame(
    level = "Protein1",
    row.names = "Protein1"
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = coldata,
    rowData = rowdata
  )

  # -----------------------------
  # Run function
  # -----------------------------
  res <- plot_trait_distribution(
    trait_se     = se,
    group_col    = "Group",
    group_levels = c("A", "B"),
    trait_name   = "TraitA",
    feature      = "Protein1"
  )

  # -----------------------------
  # Assertions
  # -----------------------------
  testthat::expect_type(res, "list")
  testthat::expect_true(all(c("p_hist", "p_box") %in% names(res)))

  testthat::expect_s3_class(res$p_hist, "ggplot")
  testthat::expect_s3_class(res$p_box, "ggplot")

  # sanity: underlying data used by ggplot exists
  testthat::expect_true(nrow(res$p_hist$data) > 0)
  testthat::expect_true(nrow(res$p_box$data) > 0)
})
