#' Import pGlyco3 GPSM results as a long-format table
#'
#' Read a pGlyco3 GPSM result file and convert it into a long-format
#' protein–peptide–glycan table with spectral counts per raw file.
#' Each row corresponds to a unique combination of protein, peptide,
#' glycan structure, and file.
#'
#' @details
#' This function takes a pGlyco3 GPSM file as input (typically named
#' \code{pGlycoDB-GP-FDR-Pro-Quant-Site.txt}). The following steps are performed:
#'
#' \itemize{
#'   \item Select relevant columns (\code{RawName}, \code{Proteins},
#'         \code{Peptide}, \code{PlausibleStruct}).
#'   \item Rename them to a standardized schema:
#'         \code{File}, \code{Protein}, \code{Peptide}, \code{GlycanStructure}.
#'   \item Collapse multiple protein IDs per PSM into a single
#'         pipe-separated string (e.g. \code{"P00123|P00456"}).
#'   \item Aggregate rows by \code{Protein}, \code{Peptide},
#'         \code{GlycanStructure}, and \code{File} as GPSM counts in each group.
#' }
#'
#' The output of this function is typically used as the input for
#' \code{\link{build_trait_se}}.
#'
#' @param gpsm_dir The path to the pGlyco3 GPSM output file
#' (for example, \code{"pGlycoDB-GP-FDR-Pro-Quant-Site.txt"}).
#'
#' @return
#' A data frame with one row per
#' \code{Protein}–\code{Peptide}–\code{GlycanStructure}–\code{File}
#' combination and the following columns:
#'
#' \itemize{
#'   \item \code{Protein} — protein identifier(s)
#'   \item \code{Peptide} — peptide sequence
#'   \item \code{GlycanStructure} — glycan structural annotation from pGlyco3
#'   \item \code{File} — raw file name
#'   \item \code{Count} — spectral count (number of GPSMs) for this combination
#' }
#'
#' @seealso
#' \code{\link{read_decipher_gpsm}},
#' \code{\link{build_trait_se}}
#'
#' @examples
#' \dontrun{
#' gpsm_table <- read_pGlyco3_gpsm("path/to/pGlycoDB-GP-FDR-Pro-Quant-Site.txt")
#' head(gpsm_table)
#' }
#'
#' @importFrom dplyr select group_by summarise n
#' @importFrom stringr str_extract_all
#' @importFrom magrittr %>%
#'
#' @export
read_pGlyco3_gpsm <- function(gpsm_dir) {
  input <- read.delim(gpsm_dir)
  input <- input %>% select("RawName", "Proteins", "Peptide", "PlausibleStruct")
  colnames(input) <- c("File", "Protein", "Peptide", "GlycanStructure")
  input$Protein <- input$Protein %>%
    str_extract_all("[A-Z0-9]+(?=_)") %>%
    lapply(unique) %>%
    sapply(function(x) paste(x, collapse = "|"))
  gpsm_table <- input %>%
    group_by(Protein, Peptide, GlycanStructure, File) %>%
    summarise(Count = n(), .groups = "drop")
  gpsm_table
}

#' Combine Glyco-Decipher GPSM results into a long-format table
#'
#' Read multiple Glyco-Decipher GPSM files from a folder,
#' merge them into a unified protein–peptide–glycan table, and attach
#' glycan structures (WURCS 2.0).
#'
#' @details
#' This function assumes that the input folder contains one or more
#' Glyco-Decipher GPSM files, typically named with the suffix
#' \code{"_GPSM_DatabaseGlycan.txt"}. For each file,
#' GPSM records are read and collapsed by \code{Protein}, \code{Peptide},
#' \code{GlycanID}, \code{File}, and \code{Count}.
#' All files are then combined into a single table, and glycan IDs are
#' mapped to WURCS structures via the global \code{glycanDatabase} object.
#' The final table uses a standardized glycan column name
#' \code{GlycanStructure} for compatibility with downstream functions.
#'
#' @param gpsm_folder_dir
#' The path to a folder containing Glyco-Decipher GPSM files (e.g. files ending with
#' \code{"_GPSM_DatabaseGlycan.txt"}).
#'
#' @return
#' A data frame in long format with one row per
#' \code{Protein}–\code{Peptide}–\code{GlycanStructure}–\code{File}
#' combination and the following columns:
#'
#' \itemize{
#'   \item \code{Protein} — protein identifier(s)
#'   \item \code{Peptide} — peptide sequence
#'   \item \code{GlycanStructure} — glycan structural annotation (WURCS 2.0)
#'   \item \code{File} — raw file name
#'   \item \code{Count} — spectral count (number of GPSMs) for this combination
#' }
#'
#' The returned table is designed to be passed to
#' \code{\link{build_trait_se}} for glycan trait computation.
#'
#' @seealso
#' \code{\link{read_pGlyco3_gpsm}},
#' \code{\link{build_trait_se}}
#'
#' @examples
#' \dontrun{
#' gpsm_table <- read_decipher_gpsm("path/to/decipher_folder/")
#' head(gpsm_table)
#' }
#'
#' @importFrom dplyr select group_by summarise bind_rows n
#' @importFrom stringr str_extract_all str_remove
#' @importFrom magrittr %>%
#'
#' @export
read_decipher_gpsm <- function(gpsm_folder_dir) {
  gpsms <- list.files(gpsm_folder_dir)
  file_names <- gpsms %>% str_remove("_GPSM_DatabaseGlycan\\.txt$")
  gpsm_paths <- file.path(gpsm_folder_dir, gpsms)
  file_GPSM_count <- function(path) {
    psm_count <-
      read.delim(path, header = TRUE, stringsAsFactors = FALSE) %>%
      select(Protein, Peptide, GlycanID, File) %>%
      group_by(Protein, Peptide, GlycanID, File) %>%
      summarise(Count = n(), .groups = "drop")

    psm_count$Protein <- psm_count$Protein %>%
      str_extract_all("[A-Z0-9]+(?=_)") %>%
      lapply(unique) %>%
      sapply(function(x) paste(x, collapse = "|"))

    psm_count
  }
  gpsm_list <- lapply(gpsm_paths, file_GPSM_count)
  gpsm_table <- bind_rows(gpsm_list)

  # Add wurcs to glyco-decipher columns
  ind <- match(gpsm_table$GlycanID, glycanDatabase$GlycanID)
  gpsm_table$GlycanID <- glycanDatabase$StructureInformation[ind]
  # Rename the colnames
  colnames(gpsm_table)[seq_len(3)] <- c("Protein", "Peptide", "GlycanStructure")

  gpsm_table
}
