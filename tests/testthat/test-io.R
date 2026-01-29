test_that("read_pGlyco3_gpsm works with toy example", {
  path <- system.file("extdata", "pGlyco3_gpsm_toyexample.txt",
    package = "glycoTraitR"
  )
  expect_true(file.exists(path))

  gpsm <- read_pGlyco3_gpsm(path)

  # base-R version returns data.frame
  expect_s3_class(gpsm, "data.frame")
  expect_true(all(c("Protein", "Peptide", "GlycanStructure", "File", "Count") %in% colnames(gpsm)))

  expect_gt(nrow(gpsm), 0)

  expect_type(gpsm$GlycanStructure, "character")
  expect_type(gpsm$Protein, "character")
  expect_type(gpsm$File, "character")

  # aggregate() should yield integer-like count after as.integer()
  expect_true(is.integer(gpsm$Count))
  expect_true(all(gpsm$Count >= 1L))

  # Protein should be collapsed into A-Z0-9 IDs separated by '|', no underscores left
  # (allow NA just in case)
  ok_pat <- is.na(gpsm$Protein) | grepl("^[A-Z0-9]+(\\|[A-Z0-9]+)*$", gpsm$Protein)
  expect_true(all(ok_pat))
})

test_that("read_decipher_gpsm works with toy example folder", {
  folder <- system.file("extdata", "decipher_toyexample", package = "glycoTraitR")
  expect_true(dir.exists(folder))

  gpsm <- read_decipher_gpsm(folder)

  expect_s3_class(gpsm, "data.frame")
  expect_true(all(c("Protein", "Peptide", "GlycanStructure", "File", "Count") %in% colnames(gpsm)))
  expect_gt(nrow(gpsm), 0)

  expect_type(gpsm$GlycanStructure, "character")
  expect_type(gpsm$Protein, "character")
  expect_true(is.integer(gpsm$Count))

  # Protein format after parsing
  ok_pat <- is.na(gpsm$Protein) | grepl("^[A-Z0-9]+(\\|[A-Z0-9]+)*$", gpsm$Protein)
  expect_true(all(ok_pat))

  # Ensure glycanDatabase mapping actually happened (i.e., looks like a structure string, not an ID)
  # More robust: compare at least one mapped value against glycanDatabase
  data("glycanDatabase", package = "glycoTraitR", envir = environment())

  # pick one row and verify it is in glycanDatabase structures
  expect_true(any(gpsm$GlycanStructure %in% glycanDatabase$StructureInformation))
})
