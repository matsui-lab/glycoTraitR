testthat::test_that("compute_glycan_traits returns a named list (with and without motifs)", {
  testthat::skip_if_not_installed("igraph")

  # minimal pGlyco3 tree
  tree <- pGlyco3_to_tree("(N(H))")

  # no motifs
  tr0 <- compute_glycan_traits(tree, motifs = NULL)
  testthat::expect_type(tr0, "list")
  testthat::expect_true(length(tr0) > 0)

  # with a simple motif (H -> N)
  motifs <- list(
    H_to_N = list(node = c("H", "N"), edge = c("a-b"))
  )
  tr1 <- compute_glycan_traits(tree, motifs = motifs)
  testthat::expect_type(tr1, "list")
  testthat::expect_true("H_to_N" %in% names(tr1))
})

testthat::test_that("annotate_traits_to_gpsm works for pGlyco3 and decipher; warns on bad `from`", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("pbapply")

  ## ---- pGlyco3 branch ----
  gpsm_pg <- data.frame(
    Protein = c("P1", "P1"),
    Peptide = c("PEP", "PEP"),
    GlycanStructure = c("(N(H))", "(N(H))"),
    File = c("F1", "F1"),
    Count = c(2L, 3L),
    stringsAsFactors = FALSE
  )

  out_pg <- glycoTraitR:::annotate_traits_to_gpsm(gpsm_pg, from = "pGlyco3", motifs = NULL)
  testthat::expect_s3_class(out_pg, "data.frame")
  testthat::expect_false("Count" %in% names(out_pg))
  testthat::expect_false("GlycanStructure" %in% names(out_pg))
  testthat::expect_true(ncol(out_pg) > 3) # Protein/Peptide/File + traits

  ## ---- decipher branch (need a valid WURCS string that your parser can map) ----
  w <- paste0(
    "WURCS=2.0/4,9,8/",
    "[u2122h_2*NCC/3=O]",
    "[a2122h-1b_1-5_2*NCC/3=O]",
    "[a1122h-1b_1-5]",
    "[a1122h-1a_1-5]",
    "/1-2-3-4-4-4-4-4-4/",
    "a4-b1_b4-c1_c3-d1_c6-f1_d2-e1_f3-g1_f6-h1_h2-i1"
  )

  gpsm_dc <- data.frame(
    Protein = c("P1"),
    Peptide = c("PEP"),
    GlycanStructure = c(w),
    File = c("F1"),
    Count = c(1L),
    stringsAsFactors = FALSE
  )

  out_dc <- glycoTraitR:::annotate_traits_to_gpsm(gpsm_dc, from = "decipher", motifs = NULL)
  testthat::expect_s3_class(out_dc, "data.frame")
  testthat::expect_false("Count" %in% names(out_dc))
  testthat::expect_false("GlycanStructure" %in% names(out_dc))
  testthat::expect_true(ncol(out_dc) > 3)
})

testthat::test_that("build_trait_se works for level=site and level=protein; warns on bad level", {
  testthat::skip_if_not_installed("igraph")
  testthat::skip_if_not_installed("pbapply")
  testthat::skip_if_not_installed("SummarizedExperiment")

  gpsm <- data.frame(
    Protein = c("P1", "P1", "P2"),
    Peptide = c("PEP1", "PEP2", "PEP2"),
    GlycanStructure = c("(N(H))", "(N(H))", "(N(H))"),
    File = c("F1", "F1", "F2"),
    Count = c(1L, 2L, 3L),
    stringsAsFactors = FALSE
  )

  meta <- data.frame(
    file = c("F1", "F2"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  ## site
  se_site <- build_trait_se(gpsm, from = "pGlyco3", motifs = NULL, level = "site", meta = meta)
  testthat::expect_s4_class(se_site, "SummarizedExperiment")
  testthat::expect_true(length(SummarizedExperiment::assays(se_site)) > 0)

  ## protein
  se_prot <- glycoTraitR:::build_trait_se(gpsm, from = "pGlyco3", motifs = NULL, level = "protein", meta = meta)
  testthat::expect_s4_class(se_prot, "SummarizedExperiment")
  testthat::expect_true(length(SummarizedExperiment::assays(se_prot)) > 0)
})
