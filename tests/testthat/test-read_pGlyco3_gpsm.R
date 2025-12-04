test_that("read_pGlyco3_gpsm works with toy example", {
  # locate built-in toy file
  path <- system.file("extdata", "pGlyco3_gpsm_toyexample.txt",
    package = "glycoTraitR"
  )

  # ensure file exists
  expect_true(file.exists(path))

  # read the file
  gpsm <- read_pGlyco3_gpsm(path)

  # basic structure checks
  expect_s3_class(gpsm, "tbl_df")
  expect_true(all(c("Protein", "Peptide", "GlycanStructure", "File", "Count") %in% colnames(gpsm)))

  # expect at least one glycoPSM
  expect_gt(nrow(gpsm), 0)

  # example glycan structure should be parsed as character vector
  expect_type(gpsm$GlycanStructure, "character")
})
