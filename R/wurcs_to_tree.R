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
  res_raw[grepl("^\\[[A-Za-z]+1122[hmok]\\b", res_raw)] <- "H"  # Mannose

  # Galactose (exclude any hexosamine with 2* N-acetyl marker)
  res_raw[grepl("^\\[[A-Za-z]+2112[hmok]\\b", res_raw) &
            !grepl("_2\\*NCC(O)?/3=O", res_raw)] <- "H"  # Galactose

  # Fucose
  res_raw[grepl("^\\[[A-Za-z]+1221[mo]\\b", res_raw)] <- "F"  # Fucose

  # HexNAc (GlcNAc, GalNAc) – N‑acetyl on C2
  res_raw[grepl("_2\\*NCC(O)?/3=O", res_raw)] <- "N"  # HexNAc

  ## ------------------- Sialic acids -------------------------------------
  # Neu5Gc : N‑glycolyl on C5
  res_raw[grepl("_5\\*NCCO/3=O", res_raw)] <- "G"  # Neu5Gc

  # Neu5Ac : N‑acetyl on C5 (but not NCCO)
  res_raw[grepl("_5\\*NCC/3=O", res_raw) &
            !grepl("_5\\*NCCO/3=O", res_raw)] <- "A"  # Neu5Ac

  ## RES-Sequence -------
  res_idx  <- as.integer(strsplit(parts[3], "-", fixed = TRUE)[[1]])
  seq_raw  <- res_raw[res_idx]

  ## LIN --------
  lin_vec <- unlist(strsplit(parts[4], "_", fixed = TRUE))
  lin_vec <- unlist(strsplit(lin_vec, "_", fixed = TRUE))
  lin_vec <- unlist(strsplit(lin_vec, "\\|"))
  lin_vec <- grep("^[a-z][0-9?]+-[a-z][0-9?]+$", lin_vec, value = TRUE)
  lin_vec <- sub("^([a-z])[0-9?]+-([a-z])[0-9?]+$", "\\1-\\2", lin_vec, perl = TRUE)

  list(seq_raw = seq_raw, lin_vec = lin_vec)
}


