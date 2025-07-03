parser_wurcs_list <- function(w){
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

  ## 1. N ︱HexNAc：have “_2*NCC/3=O” or “_2*NCCO/3=O”
  res_raw[grepl("_2\\*NCC(O)?/3=O", res_raw)] <- "N"

  ## 2. F ︱Fucose：with 1221，ended with m / o
  res_raw[grepl("^\\[[A-Za-z]+1221[mo]\\b", res_raw)] <- "F"

  ## 3. A / G ︱Sialic acids：start with A… ，have NCCO/3=O
  is_sia <- grepl("^\\[A[a-z]*[0-9]{5}[hmok]\\b", res_raw)
  res_raw[ is_sia & !grepl("NCCO/3=O", res_raw)] <- "A"   # NeuAc
  res_raw[ is_sia &  grepl("NCCO/3=O", res_raw)] <- "G"   # NeuGc

  ## 4. H ︱Hexoses（Man / Gal）：1122… or 2112… but is not HexNAc
  res_raw[grepl("^\\[[A-Za-z]+1122[hmok]\\b", res_raw)] <- "H"
  res_raw[grepl("^\\[[A-Za-z]+2112[hmok]\\b", res_raw) &
              !grepl("_2\\*NCC(O)?/3=O", res_raw)]      <- "H"

  res_raw
  ## RES-Sequence -------
  res_idx  <- as.integer(strsplit(parts[3], "-", fixed = TRUE)[[1]])
  seq_raw  <- res_raw[res_idx]

  ## LIN --------
  lin_vec <- unlist(strsplit(parts[4], "_", fixed = TRUE))
  lin_vec <- unlist(strsplit(lin_vec, "_", fixed = TRUE))
  # lin_vec <- unlist(strsplit(lin_vec, "\\|"))
  lin_vec <- sub("^([a-z])[0-9?]+-([a-z])[0-9?]+.*$", "\\1-\\2", lin_vec, perl = TRUE)

  list(seq_raw = seq_raw, lin_vec = lin_vec)
}

build_bracket <- function(lst) {
  seq_raw <- lst$seq_raw
  lin_vec <- lst$lin_vec

  n_nodes  <- length(seq_raw)
  labels   <- letters[seq_len(n_nodes)]
  mono_map <- setNames(seq_raw, labels)

  adj <- lapply(labels, function(x) character())
  names(adj) <- labels
  for (edge in lin_vec) {
    parts <- strsplit(edge, "-", fixed = TRUE)[[1]]
    a <- parts[1];  b <- parts[2]
    adj[[a]] <- c(adj[[a]], b)
    adj[[b]] <- c(adj[[b]], a)
  }

  parent <- setNames(rep(NA_character_, n_nodes), labels)
  parent["a"] <- ""                     # sentinel
  queue <- "a"
  while (length(queue)) {
    node  <- queue[1]; queue <- queue[-1]
    for (nb in adj[[node]]) {
      if (is.na(parent[nb])) {
        parent[nb] <- node
        queue      <- c(queue, nb)
      }
    }
  }

  children <- lapply(labels, function(x) character())
  names(children) <- labels
  for (ch in labels[parent != "" & !is.na(parent)]) {
    par <- parent[ch]
    children[[par]] <- c(children[[par]], ch)
  }
  children <- lapply(children, sort)

  recurse <- function(node) {
    expr <- mono_map[[node]]
    for (ch in children[[node]]) {
      expr <- paste0(expr, "(", recurse(ch), ")")
    }
    expr
  }
  paste0("(", recurse("a"), ")")
}

p <- length(glycanDatabase$StructureInformation)
r <- sapply(1:p, function(i) {
  # glycanDatabase[i,] %>% as.matrix()
  if(i %% 500 == 0) print(i);
  w <- glycanDatabase$StructureInformation[i];
  list <- parser_wurcs_list(w);
  build_bracket(list);
})

write.table(r, file = "pGlyco3_Decipher_glycanDatabase.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

