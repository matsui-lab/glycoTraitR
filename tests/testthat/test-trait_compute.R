test_that("compute_glycan_traits returns trait values", {
  tree <- list(
    node = c("N", "H"),
    edge = "a-b"
  )

  res <- compute_glycan_traits(
    tree,
    motifs = NULL
  )

  expect_type(res, "list")
  expect_gt(length(res), 0L)
})


test_that("compute_glycan_traits includes user-defined motifs", {
  tree <- list(
    node = c("N", "H"),
    edge = "a-b"
  )

  motifs <- list(
    N_to_H = list(
      node = c("N", "H"),
      edge = "a-b"
    )
  )

  res <- compute_glycan_traits(
    tree,
    motifs = motifs
  )

  expect_type(res, "list")
  expect_true("N_to_H" %in% names(res))
})
