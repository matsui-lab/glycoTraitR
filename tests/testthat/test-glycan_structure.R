test_that("build_glycan_igraph builds a directed igraph with expected attributes", {
  testthat::skip_if_not_installed("igraph")

  tree <- list(
    node = c("N", "H", "F"),
    edge = c("a-b", "a-c")
  )

  g <- build_glycan_igraph(tree)

  testthat::expect_s3_class(g, "igraph")
  testthat::expect_true(igraph::is_directed(g))

  # vertex attributes
  testthat::expect_identical(igraph::V(g)$name, c("a", "b", "c"))
  testthat::expect_identical(igraph::V(g)$residue, c("N", "H", "F"))
  testthat::expect_identical(igraph::V(g)$type, c("N", "H", "F"))
  testthat::expect_identical(igraph::V(g)$is_root, c(TRUE, FALSE, FALSE))

  # edges should be oriented away from root "a"
  ed <- igraph::as_edgelist(g, names = TRUE)
  testthat::expect_true(any(ed[, 1] == "a" & ed[, 2] == "b"))
  testthat::expect_true(any(ed[, 1] == "a" & ed[, 2] == "c"))
})

test_that("count_residues returns correct composition counts", {
  tree <- list(
    node = c("N", "H", "F", "H", "A", "G"),
    edge = c("a-b", "a-c", "b-d", "b-e", "e-f")
  )

  res <- count_residues(tree)

  testthat::expect_named(
    res,
    c("GlycanSize", "Hexose", "HexNAc", "Neu5Ac", "Neu5Gc", "Fucose")
  )

  testthat::expect_equal(unname(res["GlycanSize"]), 6)
  testthat::expect_equal(unname(res["Hexose"]), 2)
  testthat::expect_equal(unname(res["HexNAc"]), 1)
  testthat::expect_equal(unname(res["Neu5Ac"]), 1)
  testthat::expect_equal(unname(res["Neu5Gc"]), 1)
  testthat::expect_equal(unname(res["Fucose"]), 1)
})

test_that("compute_structural_traits computes expected built-in structural traits", {
  testthat::skip_if_not_installed("igraph")

  # Design a small tree that triggers:
  # - CoreFuc = 1 and AntFuc = 1 (two F under root a)
  # - Bisect = 1 (coreman has 3 kids: 2 H + 1 N)
  # - Antennas = 2 (two branch H each has an N child)
  # - Complex = 1 (>=2 antenna roots, all N)
  # - HighMan = 0, Hybrid = 0
  tree <- list(
    node = c(
      "N", # a root
      "F", # b (core fucose under root)
      "F", # c (extra F under root -> antennary fucose)
      "H", # d coreman (first H)
      "H", # e
      "H", # f
      "N", # g (bisecting N under coreman)
      "N", # h (antenna root 1)
      "N" # i (antenna root 2)
    ),
    edge = c(
      "a-b", "a-c", "a-d", # root children
      "d-e", "d-f", "d-g", # coreman children: H, H, N  -> bisecting
      "e-h", "f-i" # two antennas: H->N, H->N
    )
  )

  tr <- compute_structural_traits(tree)

  testthat::expect_type(tr, "double")
  testthat::expect_named(tr, c("Antennas", "Bisect", "Complex", "HighMan", "Hybrid", "CoreFuc", "AntFuc"))

  testthat::expect_identical(unname(tr["CoreFuc"]), 1)
  testthat::expect_identical(unname(tr["AntFuc"]), 1)
  testthat::expect_identical(unname(tr["Bisect"]), 1)
  testthat::expect_identical(unname(tr["Antennas"]), 2)
  testthat::expect_identical(unname(tr["Complex"]), 1)
  testthat::expect_identical(unname(tr["HighMan"]), 0)
  testthat::expect_identical(unname(tr["Hybrid"]), 0)
})

test_that("compute_userdefined_traits returns empty when motifs is NULL", {
  testthat::skip_if_not_installed("igraph")

  tree <- list(
    node = c("N", "H"),
    edge = c("a-b")
  )

  ud <- compute_userdefined_traits(tree, motifs = NULL)
  testthat::expect_true(is.null(ud))
  testthat::expect_length(ud, 0)
})

test_that("compute_userdefined_traits counts motif occurrences via subgraph isomorphism", {
  testthat::skip_if_not_installed("igraph")

  # same structural tree as previous structural-traits test
  tree <- list(
    node = c("N", "F", "F", "H", "H", "H", "N", "N", "N"),
    edge = c("a-b", "a-c", "a-d", "d-e", "d-f", "d-g", "e-h", "f-i")
  )

  # motif: H -> N (should match twice: e->h and f->i)
  motif_HN <- list(
    node = c("H", "N"),
    edge = c("a-b")
  )

  motifs <- list(H_to_N = motif_HN)

  ud <- compute_userdefined_traits(tree, motifs)

  testthat::expect_named(ud, "H_to_N")
  testthat::expect_equal(unname(ud["H_to_N"]), 3)
})

test_that("compute_structural_traits covers coreman-NA branch (no H present)", {
  testthat::skip_if_not_installed("igraph")

  # No "H" in node -> coreman becomes NA -> else branch should run
  tree <- list(
    node = c("N", "F", "N", "A", "G"),
    edge = c("a-b", "a-c", "c-d", "c-e")
  )

  tr <- compute_structural_traits(tree)

  testthat::expect_named(
    tr,
    c("Antennas", "Bisect", "Complex", "HighMan", "Hybrid", "CoreFuc", "AntFuc")
  )

  # else branch sets these to 0
  testthat::expect_identical(unname(tr["Antennas"]), 0)
  testthat::expect_identical(unname(tr["Bisect"]), 0)
  testthat::expect_identical(unname(tr["Complex"]), 0)
  testthat::expect_identical(unname(tr["HighMan"]), 0)
  testthat::expect_identical(unname(tr["Hybrid"]), 0)

  # CoreFuc/AntFuc still computed from root neighbors
  testthat::expect_true(unname(tr["CoreFuc"]) %in% c(0, 1))
  testthat::expect_true(unname(tr["AntFuc"]) %in% c(0, 1))
})
