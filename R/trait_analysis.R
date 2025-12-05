#' Differential analysis of glycan traits between experimental groups
#'
#' Perform group-wise statistical testing on glycan trait matrices stored in a
#' SummarizedExperiment created by \code{\link{build_trait_se}}. For each glycan
#' trait and each site/protein row, the function compares trait intensities
#' across user-specified experimental groups using Welch’s t-test and
#' Levene’s variance test.
#'
#' @details
#' Each assay in `trait_se` represents a glycan trait matrix.
#' The rows as glycopeptides (site-level) or proteins (protein-level).
#' The columns are GPSM count found in samples.
#' For each trait × feature combination,
#' Extract PSM-level trait intensities for samples belonging to the specified `group_levels`.
#' Exclude traits where either group has fewer than `min_psm` GPSMs.
#' Exclude all-zero traits (boolean-like traits)
#' Run Welch two-sample t-test (`t.test`) and
#' Levene's variance test (`car::leveneTest`, median centered).
#' A result is returned only if either test shows p < 0.05.
#'
#' @param trait_se
#' A SummarizedExperiment containing trait matrices (one assay per trait),
#' typically returned by \code{\link{build_trait_se}}.
#'
#' @param group_col
#' The column name in `colData(trait_se)` defining sample group membership.
#'
#' @param group_levels
#' Character vector specifying which values of `group_col` to compare (e.g., `c("Control","Treatment")`).
#'
#' @param min_psm
#' Minimum required PSM count per group for statistical testing.
#' Default = 20.
#'
#' @return
#' A data frame of significant trait–site (or trait–protein) comparisons with:
#'  \itemize{
#'     \item trait  — glycan trait name
#'     \item level  — site/protein identifier
#'     \item l_pval — Levene test p-value
#'     \item f_val  — Levene test F statistic
#'     \item t_pval — Welch t-test p-value
#'     \item t_val  — t-statistic
#'   }
#' Rows correspond only to significant comparisons (p < 0.05) for `l_pval` or `t_pval`.
#'
#'
#' @examples
#' # Load toy pGlyco3 GPSM data included with the package
#' path <- system.file("extdata", "pGlyco3_gpsm_toyexample.txt",
#'   package = "glycoTraitR"
#' )
#' gpsm_toyexample <- read_pGlyco3_gpsm(path)
#'
#' # Load accompanying toy metadata
#' data("meta_toyexample")
#'
#' # Build glycan trait SummarizedExperiment at the protein level
#' trait_se <- build_trait_se(
#'   gpsm = gpsm_toyexample,
#'   from = "pGlyco3",
#'   motifs = NULL,
#'   level = "protein",
#'   meta = meta_toyexample
#' )
#'
#' # Identify glycan traits significantly changed between groups
#' changed_traits <- analyze_trait_changes(
#'   trait_se     = trait_se,
#'   group_col    = "Diagnosis",
#'   group_levels = c("Normal", "Symptomatic"),
#'   min_psm      = 20
#' )
#' changed_traits
#'
#' @export
analyze_trait_changes <- function(trait_se, group_col, group_levels, min_psm = 20) {
  trait_list <- assays(trait_se)
  col_data <- colData(trait_se)
  row_data <- rowData(trait_se)

  n_level <- nrow(trait_list[[1]])
  trait_len <- length(trait_list)

  if (!group_col %in% names(col_data)) {
    stop(sprintf(
      "Column '%s' was not found in colData.\nAvailable columns: %s",
      group_col, paste(names(col_data), collapse = ", ")
    ), call. = FALSE)
  }

  col_ind <- col_data[[group_col]] %in% group_levels
  cur_group <- col_data[[group_col]][col_ind]
  pairs <- expand.grid(i = seq_len(n_level), j = seq_len(trait_len))

  trait_testing <- function(k) {
    i <- pairs$i[k]
    j <- pairs$j[k]

    trait_vec <- as.numeric(trait_list[[j]][i, col_ind])
    na_ind <- is.na(trait_vec)
    value_x <- trait_vec[!na_ind]
    meta_x <- cur_group[!na_ind]

    # skip if one of the group has less psms than threshold
    if (any(table(meta_x)[group_levels] < min_psm)) {
      return(NULL)
    }
    # skip all-zero vectors (for bool type traits)
    if (all(value_x == 0)) {
      return(NULL)
    }

    t_res <- t.test(value_x ~ meta_x, var.equal = FALSE)
    l_res <- car::leveneTest(value_x ~ meta_x, center = median)

    l_p <- l_res$`Pr(>F)`[1]
    t_p <- t_res$p.value

    if (l_p < 0.05 || t_p < 0.05) {
      return(data.frame(
        trait  = names(trait_list)[j],
        level  = row_data$level[i],
        l_pval = l_p,
        f_val  = l_res$`F value`[1],
        t_pval = t_p,
        t_val  = t_res$statistic
      ))
    } else {
      return(NULL)
    }
  }
  res_list <- pbapply::pblapply(seq_len(nrow(pairs)), trait_testing)
  do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])
}
