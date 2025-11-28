get_trait <- function(list, user_defined_traits = NULL) {

  cnt_trait <- cnt_table(list)
  struct_trait <- struct_table(list, user_defined_traits)

  traits <- c(cnt_trait, struct_trait)
  as.list(traits)
}

add_traits_to_gpsm <- function(gpsm, from, user_defined_traits = NULL) {

  # Required GPSM columns
  col_sel <- c("Protein", "Peptide", "GlycanStructure")

  # Split key GPSM columns and other metadata columns
  glycan  <- unique(gpsm$GlycanStructure)

  # Select appropriate glycan-to-tree converter
  if (from == "decipher") {
    get_trait_vector <- function(glycan) {
      glycan_list <- wurcs_to_tree(glycan)
      res_list <- get_trait(glycan_list, user_defined_traits)
      unlist(res_list)
    }
  } else if (from == "pGlyco3") {
    get_trait_vector <- function(glycan) {
      glycan_list <- pGlyco3_to_tree(glycan)
      res_list <- get_trait(glycan_list, user_defined_traits)
      unlist(res_list)
    }
  } else {
    warning("`from` must be either 'decipher' or 'pGlyco3'. No traits were computed.")
  }

  # Compute traits for each glycan structure
  traits_list <- pbapply::pblapply(glycan, get_trait_vector)
  traits_columns <- bind_rows(traits_list)

  # add trait to psm matrix
  ind <- match(gpsm$GlycanStructure, glycan)
  traits_add <- gpsm$count * traits_columns[ind,]
  # combine the both
  gpsm <- select(gpsm, -count, -GlycanStructure)
  gpsm <- cbind(gpsm, traits_add)

  gpsm
}

#' Construct a SummarizedExperiment Object Containing Glycan Trait Matrices
#'
#' @title
#' Generate peptide- or protein-level glycan trait matrices as a SummarizedExperiment.
#'
#' @description
#' This function aggregates glycan trait annotations from a GPSM table and
#' constructs a `SummarizedExperiment` object containing trait-by-PSM matrices.
#' Glycan traits are first computed at the GPSM level,
#' and are then summarized at either the glycosylation *site* level (peptide-based)
#' or *protein* level. Each trait is represented as a matrix whose rows correspond
#' to sites or proteins, and whose columns correspond to PSM identifiers.
#'
#' @details
#' The input `gpsm` must have glycopeptide-spectrum-match (GPSM) records containing:
#'
#' * **Protein** — protein identifiers
#' * **Peptide** — peptide sequence (used as glycosylation site proxy)
#' * **File** — sample or raw file name for linking PSMs to metadata
#' * Any additional GPSM-level metadata
#'
#' Glycan traits are first appended using [`add_traits_to_gpsm()`], which parses
#' glycan structures according to the search engine specification (`from`) and
#' computes user-defined or default structural traits. Each PSM is assigned a
#' unique `psm_id` to ensure correct construction of the assay matrices.
#'
#' ## Trait Matrix Generation
#' For each glycan trait column:
#'
#' 1. The GPSM table is reshaped into a long format indexed by `psm_id`.
#' 2. A wide-format matrix is created where **rows correspond to either:**
#'    * peptides (`level = "site"`), or
#'    * proteins (`level = "protein"`)
#' 3. Columns correspond to unique PSM identifiers.
#' 4. Values represent trait intensities for each PSM.
#'
#' ## Output Structure
#' The returned `SummarizedExperiment` contains:
#'
#' * **assays** — a list of matrices, one per glycan trait
#' * **rowData** — a data frame containing peptide or protein names
#' * **colData** — sample-level metadata aligned to PSMs
#'
#' This format allows downstream operations such as differential analysis,
#' tensor conversion, visualization, and multi-level trait modeling.
#'
#' @param gpsm
#' A data frame or tibble containing GPSM-level glycopeptide information.
#' Must include `Protein`, `Peptide`, `File`, and columns required by
#' `add_traits_to_gpsm()`.
#'
#' @param from
#' A character string specifying the glycan annotation format used in the GPSM
#' input. Valid values:
#'
#' * `"decipher"` — WURCS-formatted glycan structures
#' * `"pGlyco3"` — pGlyco-format glycan structures
#'
#' @param user_defined_traits
#' A named list specifying user-defined glycan structural motifs for trait
#' detection. See `add_traits_to_gpsm()` for the structure and examples.
#'
#' @param level
#' A character string indicating the summarization level:
#'
#' * `"site"` — glycan traits summarized at peptide level
#' * `"protein"` — glycan traits summarized at protein level
#'
#' Determines the row dimension of each assay matrix.
#'
#' @param meta
#' A data frame describing sample-level metadata. Must contain a `file` column
#' whose entries match the `File` column of `gpsm`. Additional metadata columns
#' are preserved and attached to the resulting `colData`.
#'
#' @return
#' A `SummarizedExperiment` object containing:
#'
#' * **assays** — named list of trait matrices (`trait × PSM`)
#' * **rowData** — data frame containing peptides or proteins used as row labels
#' * **colData** — metadata aligned to PSM identifiers
#'
#' Row and column names follow the structure required for downstream tensor
#' construction and trait modeling.
#'
#' @seealso
#' * [`add_traits_to_gpsm()`] — compute GPSM-level glycan traits
#' * `SummarizedExperiment::SummarizedExperiment`
#'
#' @examples
#' \dontrun{
#' # Example GPSM table
#' gpsm <- data.frame(
#'   Protein = c("P1", "P2"),
#'   Peptide = c("PEPTIDE1", "PEPTIDE2"),
#'   GlycanStructure = c("WURCS=...", "WURCS=..."),
#'   count = c(3, 1),
#'   File = c("rawA", "rawB")
#' )
#'
#' # Example metadata
#' meta <- data.frame(
#'   file = c("rawA", "rawB"),
#'   condition = c("A", "B")
#' )
#'
#' # Example user-defined traits
#' user_traits <- list(
#'   MotifA = list(node = c("H","N"), edge = "a-b")
#' )
#'
#' # Construct peptide-level trait SE
#' se <- obtain_trait_se(
#'   gpsm,
#'   from = "decipher",
#'   user_defined_traits = user_traits,
#'   level = "site",
#'   meta = meta
#' )
#'
#' se
#' }
#'
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom pbapply pblapply
#' @import SummarizedExperiment
#'
#' @export
obtain_trait_se <- function (gpsm,
                             from,
                             user_defined_traits,
                             level,
                             meta) {
  message("adding traits to gpsm matrix")
  gpsm_mat <- add_traits_to_gpsm(gpsm, from, user_defined_traits)

  if(level == "site") {
    sel <- "Peptide"
  } else if(level == "protein") {
    sel <- "Protein"
  } else {
    warning("`level` must be either 'site' or 'protein'.")
  }

  gpsm_mat$psm_id <- paste0("psm", 1:nrow(gpsm_mat))

  traits <- setdiff(colnames(gpsm_mat), c("Protein", "Peptide", "File", "psm_id"))
  traits_len <- length(traits)
  get_se_level_trait <- function(i) {
    out <- select(gpsm_mat, sel, "psm_id", traits[i])
    out <- pivot_wider(data = out, names_from = psm_id, values_from = traits[i])
    out <- as.data.frame(out)
    row.names(out) <- out[[1]]
    out[,-1]
  }
  message(paste("get", level, "trait matrix list"))
  trait_mat_list <- pbapply::pblapply(1:traits_len, get_se_level_trait)
  names(trait_mat_list) <- traits

  # Summarized Experiments
  ind <- match(gpsm_mat$File, meta$file)
  coldata <- meta[ind, ]
  coldata$psm_id <- gpsm_mat$psm_id

  rowdata <- data.frame(
    level = rownames(trait_mat_list[[1]])
  )

  se <- SummarizedExperiment(
    assays=trait_mat_list,
    rowData = rowdata,
    colData = coldata)
  se
}

