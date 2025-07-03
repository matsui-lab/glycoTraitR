find_likely_glycan <- function(gp_df) {
  gp_df %>%
    group_by(GlycanID) %>%
    summarise(mean_score = mean(GlycanScore, na.rm = TRUE),
              .groups = "drop") %>%
    slice_max(mean_score, n = 1, with_ties = FALSE) %>%
    pull(GlycanID)
}


find_glycanID <- function(gp_quant, gp_ident) {

  gp_quant <- gp_quant %>%
    mutate(
      pep_key = str_remove(Peptide, " .*$"),
      mod_key = if_else(str_detect(Peptide, " "), str_extract(Peptide, "(?<= )[^ ].*"), "" ),
      gly_key = str_remove(Glycan,  "\\+.*$")
    )

  gp_ident <- gp_ident %>%
    mutate(
      pep_key = str_remove(Peptide, " .*$"),
      mod_key = gp_ident$Modification,
      gly_key = str_remove(GlycanComposition, "\\+.*$")
    )

  lookup <- gp_ident %>% group_by(pep_key, mod_key, gly_key) %>%
    reframe(glycanID = find_likely_glycan(pick(GlycanScore, GlycanID)),
              .groups  = "drop")

  gp_quant %>% left_join(lookup, by = c("pep_key", "mod_key", "gly_key")) %>% pull(glycanID)
}

