testthat::test_that("build_glycan_igraph builds a directed igraph with expected attributes",
                    {
                      testthat::skip_if_not_installed("igraph")

                      tree <- list(node = c("N", "H", "F"),
                                   edge = c("a-b", "a-c"))

                      g <- build_glycan_igraph(tree)

                      testthat::expect_s3_class(g, "igraph")
                      testthat::expect_true(igraph::is_directed(g))

                      testthat::expect_identical(igraph::V(g)$name, c("a", "b", "c"))

                      testthat::expect_identical(igraph::V(g)$residue, c("N", "H", "F"))

                      testthat::expect_identical(igraph::V(g)$type, c("N", "H", "F"))

                      testthat::expect_identical(igraph::V(g)$is_root, c(TRUE, FALSE, FALSE))

                      ed <- igraph::as_edgelist(g, names = TRUE)

                      testthat::expect_true(any(ed[, 1] == "a" & ed[, 2] == "b"))

                      testthat::expect_true(any(ed[, 1] == "a" & ed[, 2] == "c"))
                    })


testthat::test_that("count_residues returns correct composition counts", {
  tree <- list(
    node = c("N", "H", "F", "H", "A", "G"),
    edge = c("a-b", "a-c", "b-d", "b-e", "e-f")
  )

  res <- count_residues(tree)

  testthat::expect_named(res,
                         c(
                           "GlycanSize",
                           "Hexose",
                           "HexNAc",
                           "Neu5Ac",
                           "Neu5Gc",
                           "Fucose"
                         ))

  testthat::expect_equal(unname(res["GlycanSize"]), 6)

  testthat::expect_equal(unname(res["Hexose"]), 2)

  testthat::expect_equal(unname(res["HexNAc"]), 1)

  testthat::expect_equal(unname(res["Neu5Ac"]), 1)

  testthat::expect_equal(unname(res["Neu5Gc"]), 1)

  testthat::expect_equal(unname(res["Fucose"]), 1)
})


testthat::test_that("compute_structural_traits computes expected built-in traits",
                    {
                      testthat::skip_if_not_installed("igraph")

                      tree <- list(
                        node = c("N", "F", "F", "H", "H", "H", "N", "N", "N"),
                        edge = c("a-b", "a-c", "a-d", "d-e", "d-f", "d-g", "e-h", "f-i")
                      )

                      tr <- compute_structural_traits(tree)

                      testthat::expect_named(
                        tr,
                        c(
                          "Antennas",
                          "IsBisecting",
                          "IsComplex",
                          "IsOligomannose",
                          "IsHybrid",
                          "IsC_Fucosed",
                          "IsA_Fucosed"
                        )
                      )

                      testthat::expect_equal(unname(tr["Antennas"]), 2)

                      testthat::expect_equal(unname(tr["IsBisecting"]), 1)

                      testthat::expect_equal(unname(tr["IsComplex"]), 1)

                      testthat::expect_equal(unname(tr["IsOligomannose"]), 0)

                      testthat::expect_equal(unname(tr["IsHybrid"]), 0)

                      testthat::expect_equal(unname(tr["IsC_Fucosed"]), 1)

                      testthat::expect_equal(unname(tr["IsA_Fucosed"]), 1)
                    })


testthat::test_that("compute_userdefined_traits returns NULL when motifs is NULL", {
  testthat::skip_if_not_installed("igraph")

  tree <- list(node = c("N", "H"), edge = "a-b")

  ud <- compute_userdefined_traits(tree, motifs = NULL)

  testthat::expect_null(ud)
})


testthat::test_that("compute_userdefined_traits counts motif occurrences", {
  testthat::skip_if_not_installed("igraph")

  tree <- list(
    node = c("N", "F", "F", "H", "H", "H", "N", "N", "N"),
    edge = c("a-b", "a-c", "a-d", "d-e", "d-f", "d-g", "e-h", "f-i")
  )

  motif_hn <- list(node = c("H", "N"), edge = "a-b")

  motifs <- list(H_to_N = motif_hn)

  ud <- compute_userdefined_traits(tree, motifs)

  testthat::expect_named(ud, "H_to_N")

  testthat::expect_equal(unname(ud["H_to_N"]), 3)
})


testthat::test_that("compute_structural_traits handles glycans without mannose", {
  testthat::skip_if_not_installed("igraph")

  tree <- list(
    node = c("N", "F", "N", "A", "G"),
    edge = c("a-b", "a-c", "c-d", "c-e")
  )

  tr <- compute_structural_traits(tree)

  testthat::expect_named(
    tr,
    c(
      "Antennas",
      "IsBisecting",
      "IsComplex",
      "IsOligomannose",
      "IsHybrid",
      "IsC_Fucosed",
      "IsA_Fucosed"
    )
  )

  testthat::expect_equal(unname(tr["Antennas"]), 0)

  testthat::expect_equal(unname(tr["IsBisecting"]), 0)

  testthat::expect_equal(unname(tr["IsComplex"]), 0)

  testthat::expect_equal(unname(tr["IsOligomannose"]), 0)

  testthat::expect_equal(unname(tr["IsHybrid"]), 0)

  testthat::expect_equal(unname(tr["IsC_Fucosed"]), 1)

  testthat::expect_equal(unname(tr["IsA_Fucosed"]), 0)
})
