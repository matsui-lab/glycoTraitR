testthat::test_that("wurcs_to_tree parses a simple WURCS structure", {
  w <- paste0(
    "WURCS=2.0/2,3,2/",
    "[a2122h-1x_1-5_2*NCC/3=O]",
    "[a1122h-1b_1-5]",
    "/1-2-1/",
    "a1-b1_b1-a1"
  )

  tree <- wurcs_to_tree(w)

  testthat::expect_type(tree, "list")
  testthat::expect_named(tree, c("node", "edge"))

  testthat::expect_length(tree$node, 3L)
  testthat::expect_length(tree$edge, 2L)

  testthat::expect_true(all(tree$node %in% c("N", "H", "A", "G", "F")))

  testthat::expect_true(all(grepl("^[A-Za-z]-[A-Za-z]$", tree$edge)))
})


testthat::test_that("pGlyco3_to_tree parses nodes and parent-child edges", {
  expr <- "(N(H(H)))"

  tree <- pGlyco3_to_tree(expr)

  testthat::expect_type(tree, "list")
  testthat::expect_named(tree, c("node", "edge"))

  testthat::expect_equal(tree$node, c("N", "H", "H"))

  testthat::expect_equal(tree$edge, c("a-b", "b-c"))
})
