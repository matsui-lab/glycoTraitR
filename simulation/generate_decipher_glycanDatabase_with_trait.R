library(readr)
glycanDatabase <- read_csv("data/decipher_related/glycanDatabase.csv", show_col_types = FALSE)
# foucs on the common N-glycan
glycanDatabase <- glycanDatabase %>% filter(Mass > 700)
wurcs <- glycanDatabase$StructureInformation

# add trait
library(pbapply)
res <- pblapply(wurcs, get_trait_from_wurcs)
res_rbind <- do.call(rbind, res)

# obtain the glytoucanID + trait
res_rbind <- cbind(res_rbind, GlycanID = glycanDatabase$GlycanID)
res_rbind <- as.data.frame(res_rbind)

res_mcp_decipher <- readRDS("~/glycoTraitR/simulation/result_mcp_decipher.rds")
gp_quant <- res_mcp_decipher$quant_matrix
meta <- res_mcp_decipher$meta
ind <- match(gp_quant$glycan_id, res_rbind$GlycanID)
gp_quant <- cbind(gp_quant, res_rbind[ind, ])

# quantification matrix
mat <- gp_quant[, 4:37]
mat[is.na(mat)] <- 0
mat_num <- as.matrix(mat)
heatmap(mat_num)
