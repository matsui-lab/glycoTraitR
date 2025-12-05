#' Toy metadata for glycoTraitR examples
#'
#' A toy example metadata table used in vignettes, examples, and tests.
#'
#' @format A data frame with 34 rows and 2 variables:
#' \describe{
#'   \item{Diagnosis}{Clinical diagnosis (factor)}
#'   \item{Sample number}{Sample ID (numeric)}
#'   \item{file}{File Name (character)}
#' }
#'
#' @usage data(meta_toyexample)
#' @examples
#' data(meta_toyexample)
#' head(meta_toyexample)
"meta_toyexample"

#' Glycan annotation reference database
#'
#' A curated reference table mapping glycan IDs to their structures,
#' used internally by glycoTraitR and also available to users.
#'
#' @format A data frame with N rows and M columns.
#'
#' @usage data(glycanDatabase)
#'
#' @examples
#' data(glycanDatabase)
#' head(glycanDatabase)
"glycanDatabase"
