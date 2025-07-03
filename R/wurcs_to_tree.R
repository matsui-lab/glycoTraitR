# Parse a WURCS string into its main components
parser_wurcs <- function(w){
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
  res_raw[grepl("^\\[[A-Za-z]+1122[hmok]\\b",   res_raw)] <- "H" # Mannose
  res_raw[grepl("^\\[[A-Za-z]+2112[hmok]\\b",   res_raw) &
            !grepl("_2\\*NCC(O)?/3=O",          res_raw)] <- "H" # Galactose
  res_raw[grepl("^\\[[A-Za-z]+1221[mo]\\b",     res_raw)] <- "F" # Fucose
  res_raw[grepl("_2\\*NCC(O)?/3=O",             res_raw)] <- "N" # HexNAc
  res_raw[grepl("^\\[A[a-z]*[0-9]{5}[hmok]\\b", res_raw)] <- "A" # NeuAc
  res_raw[grepl("^\\[A[a-z]*[0-9]{5}[hmok]\\b", res_raw)] <- "G" # NeuGc

  ## RES-Sequence -------
  res_idx  <- as.integer(strsplit(parts[3], "-", fixed = TRUE)[[1]])
  seq_raw  <- res_raw[res_idx]

  ## LIN --------
  lin_vec <- unlist(strsplit(parts[4], "_", fixed = TRUE))
  lin_vec <- unlist(strsplit(lin_vec, "_", fixed = TRUE))
  lin_vec <- unlist(strsplit(lin_vec, "\\|"))
  lin_vec <- grep("^[a-z][0-9?]+-[a-z][0-9?]+$", lin_vec, value = TRUE)

  list(seq_raw = seq_raw, lin_vec = lin_vec)
}
