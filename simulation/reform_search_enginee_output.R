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

base_dir <- "./data/decipher_related"
files <- list.files(paste0(base_dir,"/GPSM"))
sample_names <- str_remove(files, "_GPSM_DatabaseGlycan\\.txt$")

get_file_paths <- function(sample_name) {
  list(
    quant = file.path(base_dir, "decipher_mcp", paste0(sample_name, "_GlycoPeptideQuantificationArea.txt")),
    ident = file.path(base_dir, "GPSM", paste0(sample_name, "_GPSM_DatabaseGlycan.txt")),
    sample = sample_name
  )
}

file_list <- lapply(sample_names, get_file_paths)

all_data <- lapply(file_list, function(files) {
  gp_quant <- read.delim(files$quant, header = TRUE, stringsAsFactors = FALSE)
  gp_ident <- read.delim(files$ident, header = TRUE, stringsAsFactors = FALSE)

  gp_quant$glycan_id <- find_glycanID(gp_quant, gp_ident)

  gp_quant %>%
    select(Peptide, Site, glycan_id, Area)
})

library(tidyr)
quant_matrix <- reduce(all_data, full_join, by = c("Peptide", "Site", "glycan_id"))
file_names <- lapply(file_list, function(x) x$sample) %>% unlist
colnames(quant_matrix)[4:37] <- file_names

sample_id <- str_extract_all(file_names, "(?<=AD)[0-9]+") %>% as.numeric()

library(readxl)
meta <- read_excel("data/meta_table.xlsx")
ind <- match(sample_id, meta$`Sample number`)
meta <- meta[ind, ]

res_list <- list(
  quant_matrix = quant_matrix,
  meta = meta
)

saveRDS(res_list, file = "./simulation/result_mcp_decipher.rds")



