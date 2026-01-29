test_that("wurcs_to_tree parses nodes and edges correctly", {

  w <- paste0(
    "WURCS=2.0/2,3,2/",
    "[a2122h-1x_1-5_2*NCC/3=O]",
    "[a1122h-1b_1-5]",
    "/1-2-1/",
    "a1-b1_b1-a1"
  )

  tree <- wurcs_to_tree(w)

  ## basic structure
  expect_type(tree, "list")
  expect_true(all(c("node", "edge") %in% names(tree)))

  ## node logic
  expect_type(tree$node, "character")
  expect_gt(length(tree$node), 0)

  ## mapping via WURCS_RES_MAP happened
  expect_true(all(tree$node %in% unname(WURCS_RES_MAP)))

  ## edge parsing
  expect_type(tree$edge, "character")
  expect_true(all(grepl("^[A-Za-z]-[A-Za-z]$", tree$edge)))
})


test_that("pGlyco3_to_tree parses parentheses and edges correctly", {

  expr <- "(N(H(H)))"
  tree <- pGlyco3_to_tree(expr)

  ## basic structure
  expect_type(tree, "list")
  expect_true(all(c("node", "edge") %in% names(tree)))

  ## nodes are letters
  expect_type(tree$node, "character")
  expect_true(all(tree$node %in% c("N", "H")))

  ## edges exist and have correct format
  expect_type(tree$edge, "character")
  expect_true(all(grepl("^[a-zA-Z]-[a-zA-Z]$", tree$edge)))

  ## root has children (covers parent != NULL branch)
  expect_gt(length(tree$edge), 0)
})

