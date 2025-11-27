obtain_changed_trait_table <- function(trait_se, group_list, psm_thres = 20) {

  trait_list <- assays(trait_se)
  col_data   <- colData(trait_se)
  row_data   <- rowData(trait_se)

  n_level   <- nrow(trait_list[[1]])
  trait_len <- length(trait_list)

  # group columns
  col_sel   <- names(group_list)
  group     <- group_list[[col_sel]]
  if(! col_sel %in% names(col_data)) {
    stop(
      sprintf(
        "Column '%s' specified in group_list was not found in colData. ",
        col_sel
      ),
      "Available columns are: ",
      paste(names(col_data), collapse = ", "),
      call. = FALSE
    )
  }
  col_ind   <- col_data[[col_sel]] %in% group
  cur_group <- col_data[[col_sel]][col_ind]

  pairs <- expand.grid(i = seq_len(n_level), j = seq_len(trait_len))

  trait_testing <- function(k) {
    i <- pairs$i[k]
    j <- pairs$j[k]

    trait_vec <- as.numeric(trait_list[[j]][i, col_ind])
    na_ind <- is.na(trait_vec)
    value_x <- trait_vec[!na_ind]
    meta_x  <- cur_group[!na_ind]

    # skip if one of the group has less psms than threshold
    if (any(table(meta_x)[group] < psm_thres))
      return(NULL)
    # skip all-zero vectors (for bool type traits)
    if (all(value_x == 0))
      return(NULL)

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
    }
    return(NULL)
  }

  res_list <- pbapply::pblapply(seq_len(nrow(pairs)), trait_testing)
  res_df <- do.call(rbind, res_list[!sapply(res_list, is.null)])
  res_df
}
