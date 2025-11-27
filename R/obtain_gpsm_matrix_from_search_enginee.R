#' Extract GPSM Matrix from pGlyco3 Output
#'
#' This function reads a pGlyco3 GPSM result file and converts it into a
#' protein–peptide–glycan matrix. It extracts protein, peptide and GPSM counts,
#' and reshapes the data into a wide matrix where each file is a column.
#'
#' @param gpsm_dir Path to the pGlyco3 GPSM output (pGlycoDB-GP-FDR-Pro-Quant-Site.txt) file.
#'
#' @return A wide-format data frame where:
#' \itemize{
#'   \item \code{Protein} — Protein name
#'   \item \code{Peptide} — Peptide sequences
#'   \item \code{GlycanStructure} — Glycan structural annotations from pGlyco3
#'   \item Sample columns — One column per raw file, containing spectral counts
#' }
#'
#' @examples
#' \dontrun{
#' gpsm_matrix <- obtain_pGlyco3_gpsm("path/to/pGlyco3_GPSM.txt")
#' }
#'
#' @export
obtain_pGlyco3_gpsm <- function(gpsm_dir) {
  input <- read.delim(gpsm_dir)
  input <- input %>% select("RawName", "Proteins", "Peptide", "PlausibleStruct")
  colnames(input) <- c("File", "Protein", "Peptide", "GlycanStructure")
  input$Protein <- input$Protein %>%
    str_extract_all("[A-Z0-9]+(?=_)") %>%
    lapply(unique) %>%
    sapply(function(x) paste(x, collapse = "|"))
  gpsm_matrix <- input %>%
    group_by(Protein, Peptide, GlycanStructure, File) %>%
    summarise(count = n(), .groups = "drop")
  gpsm_matrix
}

#' Combine Glyco-Decipher GPSM Results Into a Matrix
#'
#' This function reads multiple Glyco-Decipher GPSM result files in a folder,
#' merges them by Protein–Peptide–GlycanID, converts GlycanID to
#' glycan structures (wurcs 2.0), and outputs a unified GPSM matrix where each file is a column.
#'
#' @param gpsm_folder_dir A character string. Path to a folder containing
#' Glyco-Decipher GPSM result files.
#'
#' @return A wide-format data frame with:
#' \itemize{
#'   \item \code{Protein} — Protein name
#'   \item \code{Peptide} - Peptide sequence
#'   \item \code{GlycanStructure} — Glycan structural annotations from Decipher (wurcs2.0)
#'   \item Sample columns — One column per input file with spectrum counts
#' }
#'
#' @examples
#' \dontrun{
#' gpsm_matrix <- obtain_decipher_gpsm("path/to/decipher_folder/")
#' }
#'
#' @export
obtain_decipher_gpsm <- function(gpsm_folder_dir) {
  gpsms <- list.files(gpsm_folder_dir)
  file_names <- gpsms %>% str_remove("_GPSM_DatabaseGlycan\\.txt$")
  gpsm_pathes <- file.path(gpsm_folder_dir, gpsms)
  file_GPSM_count <- function(path) {
    psm_count <-
      read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
      select(Protein, Peptide, GlycanID, File) %>%
      group_by(Protein, Peptide, GlycanID, File) %>%
      summarise(count = n(), .groups = "drop")

    psm_count$Protein <-  psm_count$Protein %>%
      str_extract_all("[A-Z0-9]+(?=_)") %>%
      lapply(unique) %>%
      sapply(function(x) paste(x, collapse = "|"))

    psm_count
  }
  gpsm_list <- lapply(gpsm_pathes, file_GPSM_count)
  gpsm_matrix <- bind_rows(gpsm_list)

  # Add wurcs to glyco-decipher columns
  ind <- match(gpsm_matrix$GlycanID, glycanDatabase$GlycanID)
  gpsm_matrix$GlycanID <- glycanDatabase$StructureInformation[ind]
  # Rename the colnames
  colnames(gpsm_matrix)[1:3] <- c("Protein", "Peptide", "GlycanStructure")

  gpsm_matrix
}
