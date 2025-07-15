res_cat <- c(
  "[a2122h-1x_1-5_2*NCC/3=O]"      = "N",
  "[a1221m-1x_1-5]"                = "F",
  "[a1122h-1x_1-5]"                = "H",
  "[a2122h-1b_1-5_2*NCC/3=O]"      = "N",
  "[a1122h-1b_1-5]"                = "H",
  "[a1122h-1a_1-5]"                = "H",
  "[a2122h-1a_1-5_2*NCC/3=O]"      = "N",
  "[u2122h_2*NCC/3=O]"             = "N",
  "[a2112h-1b_1-5]"                = "H",
  "[a1221m-1a_1-5]"                = "F",
  "[Aad21122h-2a_2-6_5*NCC/3=O]"   = "A",
  "[a2112h-1a_1-5]"                = "H",
  "[a2112h-1b_1-5_2*NCC/3=O]"      = "N",
  "[Aad21122h-2a_2-6_5*NCCO/3=O]"  = "G",
  "[a1221m-1b_1-5]"                = "F",
  "[a2112m-1a_1-5]"                = "F",
  "[o2122h_2*NCC/3=O]"             = "N",
  "[a2112h-1b_1-4]"                = "H",
  "[a4334h-1b_1-4]"                = "H",
  "[a4344h-1x_1-5_2*NCC/3=O]"      = "N",
  "[axxxxh-1x_1-5_2*NCC/3=O]"      = "N",
  "[axxxxh-1x_1-?_2*NCC/3=O]"      = "N",
  "[a2122h-1b_1-5]"                = "H",
  "[a2112h-1x_1-5]"                = "H",
  "[axxxxh-1b_1-5_2*NCC/3=O]"      = "N",
  "[a2112h-1x_1-5_2*NCC/3=O]"      = "N",
  "[axxxxh-1b_1-?_2*NCC/3=O]"      = "N",
  "[a2112h-1a_1-4]"                = "H",
  "[a2122h-1a_1-5]"                = "H",
  "[a4334h-1a_1-4]"                = "H",
  "[a4344h-1a_1-?]"                = "H",
  "[a2112m-1b_1-5]"                = "F",
  "[Aad21122h-2x_2-6_5*NCC/3=O]"   = "A",
  "[a2112h-1x_1-4]"                = "H",
  "[a4334h-1x_1-4]"                = "H",
  "[axxxxh-1x_1-5]"                = "H",
  "[a1121h-1a_1-5]"                = "H",
  "[a2122h-1x_1-5]"                = "H",
  "[a4344h-1x_1-?]"                = "H",
  "[Aad21122h-2x_2-6_5*NCCO/3=O]"  = "G",
  "[a2112h-1a_1-5_2*NCC/3=O]"      = "N",
  "[axxxxh-1x_1-?]"                = "H",
  "[a4334h-1b_1-5]"                = "H",
  "[a4344h-1b_1-5_2*NCC/3=O]"      = "N",
  "[a2122h-1b_1-?]"                = "H",
  "[a2122h-1b_1-?_2*NCC/3=O]"      = "N"
)

# Parse a WURCS string into its main components
wurcs_to_tree <- function(w){
  core <- sub("^WURCS=[^/]+/", "", w)

  parts <- character()
  cur <- ""
  depth <- 0L
  for(ch in strsplit(core, "")[[1]]){
    if (ch == "[") depth <- depth + 1L
    if (ch == "]") depth <- depth - 1L
    if (ch == "/" && depth == 0L){
      parts <- c(parts, cur)
      cur <- ""
    } else {
      cur <- paste0(cur, ch)
    }
  }
  parts <- c(parts, cur)

  ## UniqueRES --------
  res_raw  <- regmatches(parts[2],
                         gregexpr("\\[[^]]+\\]", parts[2], perl = TRUE))[[1]]
  res_raw <- res_cat[res_raw]

  ## RES-Sequence -------
  res_idx  <- as.integer(strsplit(parts[3], "-", fixed = TRUE)[[1]])
  seq_raw  <- res_raw[res_idx]

  ## LIN --------
  lin_vec <- unlist(strsplit(parts[4], "_", fixed = TRUE))
  lin_vec <- unlist(strsplit(lin_vec, "_", fixed = TRUE))
  lin_vec <- unlist(strsplit(lin_vec, "\\|"))
  lin_vec <- grep("^[A-Za-z][0-9?]+-[A-Za-z][0-9?]+$", lin_vec, value = TRUE)
  lin_vec <- sub("^([A-Za-z])[0-9?]+-([A-Za-z])[0-9?]+$", "\\1-\\2", lin_vec, perl = TRUE)

  list(seq_raw = seq_raw, lin_vec = lin_vec)
}




