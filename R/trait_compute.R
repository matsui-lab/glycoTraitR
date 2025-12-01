#' Compute glycan traits from a parsed glycan tree
#'
#' Combine residue-level composition traits (see \code{\link{count_residues}}),
#' structural traits (see \code{\link{compute_structural_traits}}), and
#' user-defined motifs (see \code{\link{compute_userdefined_traits}})
#' into a unified trait vector.
#'
#' @param tree A parsed glycan tree from \code{\link{pGlyco3_to_tree}} or \code{\link{wurcs_to_tree}}.
#'
#' @param motifs Optional named list of user-defined glycan motifs.
#'
#' @return A named list of numeric trait values.
#'
#' @keywords internal
compute_glycan_traits <- function(tree, motifs) {

  cnt_trait <- count_residues(tree)
  struct_trait <- compute_structural_traits(tree)
  ud_trait <- compute_userdefined_traits(tree, motifs)

  traits <- c(cnt_trait, struct_trait, ud_trait)
  as.list(traits)
}

#' Append computed glycan traits to a GPSM table
#'
#' Parse each unique glycan structure in a GPSM table, compute its
#' residue and structural traits, and expand the resulting trait
#' columns back to the full GPSM table.
#'
#' @param gpsm A GPSM table from either
#' \code{\link{read_pGlyco3_gpsm}} or \code{\link{read_decipher_gpsm}}.
#'
#' @param from One of `decipher` or `pGlyco3`.
#'
#' @param motifs Optional list of user-defined structural motifs.
#'
#' @return A modified GPSM table with appended trait columns.
#'
#' @keywords internal
annotate_traits_to_gpsm <- function(gpsm, from, motifs) {


  # Split key GPSM columns and other metadata columns
  glycan  <- unique(gpsm$GlycanStructure)

  # Select appropriate glycan-to-tree converter
  if (from == "decipher") {
    get_trait_vector <- function(glycan) {
      glycan_list <- wurcs_to_tree(glycan)
      res_list <- compute_glycan_traits(glycan_list, motifs)
      unlist(res_list)
    }
  } else if (from == "pGlyco3") {
    get_trait_vector <- function(glycan) {
      glycan_list <- pGlyco3_to_tree(glycan)
      res_list <- compute_glycan_traits(glycan_list, motifs)
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
  traits_add <- gpsm$Count * traits_columns[ind,]
  # combine the both
  gpsm <- select(gpsm, -Count, -GlycanStructure)
  gpsm <- cbind(gpsm, traits_add)

  gpsm
}

#' Build a SummarizedExperiment of glycan trait matrices
#'
#' Convert a GPSM table into peptide- or protein-level glycan trait matrices
#' and store them in a \code{SummarizedExperiment} object. Each trait becomes an assay
#' matrix whose rows represent peptides or proteins, and whose columns
#' represent individual GPSMs.
#' This function provides a unified container for downstream analyses such as
#' differential testing and visualization.
#'
#' @param gpsm A GPSM table containing at least:
#'   `Protein`, `Peptide`, `GlycanStructure`, `File`, and `Count`.
#' @param from Character; glycan format used in the GPSM input.
#'   One of `decipher` or `pGlyco3`.
#' @param motifs Optional named list of user-defined motif structures
#'   passed to \link{compute_glycan_traits}.
#' @param level Summarization level. Either `site (peptide)` or `"protein"`.
#' @param meta Data frame of sample metadata with a column `file`
#'   matching `gpsm$File`.
#'
#' @return A `SummarizedExperiment` where each assay is a glycan-trait
#'   matrix (trait × PSM), `rowData` contains peptide/protein names,
#'   and `colData` contains metadata aligned to PSMs.
#'
#' @importFrom tidyr pivot_wider
#' @importFrom SummarizedExperiment assays colData rowData
#'
#' @export
build_trait_se <- function (gpsm, from, motifs = NULL, level, meta) {

  message("adding traits to the gpsm matrix")
  gpsm_mat <- annotate_traits_to_gpsm(gpsm, from, motifs)

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
    out <- tidyr::pivot_wider(data = out, names_from = psm_id, values_from = traits[i])
    out <- as.data.frame(out)
    row.names(out) <- out[[1]]
    out[,-1]
  }
  message(paste("generating", level, "trait matrices"))
  trait_mat_list <- pbapply::pblapply(1:traits_len, get_se_level_trait)
  names(trait_mat_list) <- traits

  # Summarized Experiments
  ind <- match(gpsm_mat$File, meta$file)
  coldata <- meta[ind, ]
  coldata$psm_id <- gpsm_mat$psm_id

  rowdata <- data.frame(
    level = rownames(trait_mat_list[[1]])
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays=trait_mat_list,
    rowData = rowdata,
    colData = coldata)
  se
}

