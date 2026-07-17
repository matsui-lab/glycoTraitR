#' Statistical testing of glycan heterogeneity scores
#'
#' Test whether trait-wise glycan heterogeneity scores differ between two
#' experimental groups using either a permutation test or Welch's two-sample
#' t-test.
#'
#' The input should be a site- or protein-level heterogeneity table generated
#' from GPSM data. Columns ending in `_h` or `_mu` are treated as score columns
#' and tested feature by feature.
#'
#' @details
#' For each trait and each site/protein feature, the observed statistic is the
#' difference in mean score between the two groups:
#'
#' \deqn{
#' T = \bar{h}_{group1} - \bar{h}_{group2}
#' }
#'
#' The statistical test is determined by `method`.
#'
#' When `method = "permutation"`, group labels are randomly permuted `B`
#' times to generate an empirical null distribution. The two-sided
#' permutation p-value is calculated as:
#'
#' \deqn{
#' p = \frac{1 + \sum_b I(|T_b| \ge |T|)}{1 + B}
#' }
#'
#' When `method = "t.test"`, Welch's two-sample t-test is applied to the
#' two groups.
#'
#' The permutation test is distribution-free and applies to both
#' integer-valued and Boolean glycan traits.
#'
#' @param h_tab A site- or protein-level score table, typically one of
#'   `hscore$h_site` or `hscore$h_protein`. Must contain a `File` column, a
#'   feature column such as `Peptide` or `Protein`, and one or more score
#'   columns ending in `_h` or `_mu`.
#'
#' @param meta A sample metadata data frame containing `File` and the grouping
#'   column specified by `group_col`.
#'
#' @param group_col Character string specifying the metadata column that defines
#'   experimental group membership.
#'
#' @param group_levels Character vector of length 2 specifying the two groups
#'   to compare, for example `c("Normal", "Symptomatic")`. The reported
#'   difference is `group_levels[1] - group_levels[2]`.
#'
#' @param feature_col Character string specifying the feature identifier column
#'   in `h_tab`. Use `"Peptide"` for site-level testing and `"Protein"` for
#'   protein-level testing.
#'
#' @param B Integer; number of random label permutations used when
#'   `method = "permutation"`. Default is 1000.
#'
#' @param min_samples Integer; minimum number of non-missing samples required
#'   in each group for a feature-trait test. Default is 3.
#'
#'
#' @param method Statistical test used for group comparison. One of
#'   `"permutation"` or `"t.test"`. Default is `"permutation"`.
#'
#' @return A data frame with one row per tested trait-feature combination and
#'   the following columns:
#' \itemize{
#'   \item `trait`: glycan trait name, with the `_h` or `_mu` suffix removed.
#'   \item `score_type`: score type, either `"h"` or `"mu"`.
#'   \item `feature`: site or protein identifier.
#'   \item `group1`, `group2`: compared group labels.
#'   \item `mean_group1`, `mean_group2`: group mean scores.
#'   \item `diff`: observed mean difference, `group1 - group2`.
#'   \item `pval`: two-sided p-value.
#'   \item `n_group1`, `n_group2`: number of non-missing samples per group.
#'   \item `method`: statistical test used.
#' }
#'
#' @keywords internal
#' @noRd
test_hscore_changes <- function(h_tab,
                                meta,
                                group_col,
                                group_levels,
                                feature_col,
                                B = 1000,
                                min_samples = 3,
                                method = c("permutation", "t.test")) {
  method <- match.arg(method)

  score_cols <- grep("_(h|mu)$", colnames(h_tab), value = TRUE)

  meta_sub <- meta[meta[[group_col]] %in% group_levels,
                   c("File", group_col), drop = FALSE]

  dat <- merge(h_tab, meta_sub, by.x = "File", by.y = "File")
  dat[[group_col]] <- factor(dat[[group_col]], levels = group_levels)

  features <- unique(dat[[feature_col]])
  pairs <- expand.grid(
    feature = features,
    trait = score_cols,
    stringsAsFactors = FALSE
  )

  one_test <- function(k) {
    feature_k <- pairs$feature[k]
    trait_k <- pairs$trait[k]

    sub <- dat[dat[[feature_col]] == feature_k, , drop = FALSE]

    value <- as.numeric(sub[[trait_k]])
    group <- sub[[group_col]]

    keep <- !is.na(value) & !is.na(group)
    value <- value[keep]
    group <- group[keep]

    n1 <- sum(group == group_levels[1])
    n2 <- sum(group == group_levels[2])
    if (n1 < min_samples || n2 < min_samples) return(NULL)

    x1 <- value[group == group_levels[1]]
    x2 <- value[group == group_levels[2]]
    mean1 <- mean(x1)
    mean2 <- mean(x2)
    obs_diff <- mean1 - mean2

    if (length(unique(c(x1, x2))) == 1) {
      return(NULL)
    }

    if (method == "permutation") {
      perm_diff <- replicate(B, {
        perm_group <- sample(group, length(group), replace = FALSE)
        mean(value[perm_group == group_levels[1]]) -
          mean(value[perm_group == group_levels[2]])
      })

      pval <- (1 + sum(abs(perm_diff) >= abs(obs_diff))) / (1 + B)

    } else if (method == "t.test") {
      pval <- tryCatch(
        stats::t.test(x1, x2)$p.value,
        error = function(e) NA_real_
      )

      if (is.na(pval)) return(NULL)
    }

    data.frame(
      trait = sub("_(h|mu)$", "", trait_k),
      score_type = sub("^.*_(h|mu)$", "\\1", trait_k),
      feature = feature_k,
      group1 = group_levels[1],
      group2 = group_levels[2],
      mean_group1 = mean1,
      mean_group2 = mean2,
      diff = obs_diff,
      pval = pval,
      n_group1 = n1,
      n_group2 = n2,
      method = method,
      stringsAsFactors = FALSE
    )
  }

  res_list <- pbapply::pblapply(seq_len(nrow(pairs)), one_test)
  res <- do.call(
    rbind,
    res_list[!vapply(res_list, is.null, logical(1))]
  )

  rownames(res) <- NULL
  res
}



#' Differential testing of site- and protein-level H-scores
#'
#' Apply statistical testing to compare site-level and protein-level glycan
#' heterogeneity scores between two experimental groups, then combine the
#' results into a single summary table.
#'
#' @param hscore A list containing `h_site` and `h_protein`, typically returned
#'   by \code{\link{compute_hscore}}.
#' @param meta A sample metadata data frame containing `File` and the grouping
#'   column specified by `group_col`.
#' @param group_col Character string specifying the metadata column defining
#'   experimental group membership.
#' @param group_levels Character vector of length 2 specifying the two groups
#'   to compare. The reported difference is `group_levels[1] - group_levels[2]`.
#' @param B Integer; number of random label permutations used when
#'   `method = "permutation"`. Default is 1000.
#' @param min_samples Integer; minimum number of non-missing samples required
#'   in each group for each feature-trait test. Default is 3.
#' @param method Statistical test used for group comparison. One of
#'   `"permutation"` or `"t.test"`. Default is `"permutation"`.
#'
#' @return A data frame combining site- and protein-level test results. The
#'   column `level` indicates whether each row is from site-level (`"site"`)
#'   or protein-level (`"protein"`) testing.
#'
#' @keywords internal
#' @noRd
test_hscore_changes_all <- function(hscore,
                                    meta,
                                    group_col,
                                    group_levels,
                                    B = 1000,
                                    min_samples = 3,
                                    method = c("permutation", "t.test")) {
  method <- match.arg(method)

  message("Test heterogeneity change at site")
  site_res <- test_hscore_changes(
    h_tab = hscore$h_site,
    meta = meta,
    group_col = group_col,
    group_levels = group_levels,
    feature_col = "Peptide",
    B = B,
    min_samples = min_samples,
    method = method
  )
  message("Test heterogeneity change at protein")
  protein_res <- test_hscore_changes(
    h_tab = hscore$h_protein,
    meta = meta,
    group_col = group_col,
    group_levels = group_levels,
    feature_col = "Protein",
    B = B,
    min_samples = min_samples,
    method = method
  )

  site_res$level <- "site"
  protein_res$level <- "protein"
  res <- rbind(site_res, protein_res)
  rownames(res) <- NULL
  res
}


#' Analyze differential glycan trait means and heterogeneity
#'
#' @description
#' A high-level wrapper that computes glycan trait abundance and heterogeneity
#' scores from GPSM data and tests for differences between two experimental
#' groups.
#'
#' The function internally:
#' \enumerate{
#'   \item Computes site- and protein-level glycan heterogeneity scores.
#'   \item Performs permutation-based or Welch's t-test group comparisons.
#' }
#'
#' @param gpsm A GPSM table imported by
#'   \code{\link{read_pGlyco3_gpsm}} or
#'   \code{\link{read_decipher_gpsm}}.
#' @param from Character string specifying the glycan representation.
#'   One of \code{"pGlyco3"} or \code{"decipher"}.
#' @param motifs Optional glycan motif definitions.
#' @param meta Sample metadata table.
#' @param group_col Column name in \code{meta} defining the comparison groups.
#' @param group_levels Length-2 character vector specifying the two groups
#'   to compare.
#' @param B Integer; number of random label permutations used when
#'   `method = "permutation"`. Default is 1000.
#' @param min_samples Minimum number of non-missing samples required in each
#'   group.
#' @param method Statistical test used for group comparison. One of
#'   `"permutation"` or `"t.test"`. Default is `"permutation"`.
#'
#' @return
#' A data frame containing differential analysis results for all glycan traits.
#'
#' @examples
#' gpsm <- readRDS(system.file("extdata", "gpsm_toyexample.rds", package = "glycoTraitR"))
#' meta <- readRDS(system.file("extdata", "meta_toyexample.rds", package = "glycoTraitR"))
#'
#' res <- analyze_hscore_changes(
#'   gpsm = gpsm,
#'   from = "pGlyco3",
#'   meta = meta,
#'   group_col = "Diagnosis",
#'   group_levels = c("Normal", "Symptomatic"),
#'   B = 5
#' )
#'
#' head(res)
#'
#' @export
analyze_hscore_changes <- function(gpsm,
                                   from,
                                   motifs = NULL,
                                   meta,
                                   group_col,
                                   group_levels,
                                   B = 1000,
                                   min_samples = 3,
                                   method = c("permutation", "t.test")) {
  method <- match.arg(method)

  message("Compute heterogeneity scores at each site and protein")
  hscore <- compute_hscore(
    gpsm = gpsm,
    from = from,
    motifs = motifs
  )

  res <- test_hscore_changes_all(
    hscore = hscore,
    meta = meta,
    group_col = group_col,
    group_levels = group_levels,
    B = B,
    min_samples = min_samples,
    method = method
  )

  res
}
