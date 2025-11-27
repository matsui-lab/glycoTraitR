library(dplyr)
library(readxl)

# Glyco-Decipher Input Files
gpsm_dir <- "./simulation/data/decipher_related/GPSM_mcp_2022_100433/"
files <- list.files(gpsm_dir)
meta_id <- files %>% stringr::str_extract("(?<=KS_AD)\\d+") %>% as.integer()

# Glyco-Decipher colData File
meta  <- read_xlsx("./simulation/data/meta_mcp_2022_100433.xlsx", col_names = TRUE)
ind <- match(meta_id, meta$`Sample number`) # No 27 and No 24 have three shots
meta <- meta[ind, ]
meta$Diagnosis <- factor(meta$Diagnosis, levels = c("Normal", "Asymptomatic", "Symptomatic"))
meta$file <- files %>% stringr::str_remove("_GPSM_DatabaseGlycan\\.txt$")

saveRDS(meta, "./simulation/meta_mcp_2022_100433.rds")
