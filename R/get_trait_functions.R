

man_cnt    <- function(list) sum(list$seq_raw == "Mannose")
gal_cnt    <- function(list) sum(list$seq_raw == "Galactose")
fuc_cnt    <- function(list) sum(list$seq_raw == "Fucose")
hexnac_cnt <- function(list) sum(list$seq_raw == "HexNAc")
sial_cnt   <- function(list) sum(list$seq_raw == "Sialic")

is_corefuc <- function(list) {
  fuc_mask <- list$seq_raw == "Fucose"
  res_letters <- letters[seq_along(list$seq_raw)]
  fuc_letters <- res_letters[fuc_mask]
  if(length(fuc_letters) == 0) return(FALSE)
  else {
    core_fuc_letters <- fuc_letters[
      sapply(fuc_letters, function(fl){
        any(grepl(paste0("^a6-",fl,"1$"), list$lin_vec)) |
          any(grepl(paste0("^", fl,"1-a6$"), list$lin_vec))
      })
    ]
    return(length(core_fuc_letters))
  }
}

sial_link_cnt <- function(list) {
  sial_mask <- list$seq_raw == "Sialic"
  res_letters <- letters[seq_along(list$seq_raw)]
  sial_letters <- res_letters[sial_mask]

  cnt23 <- cnt26 <- 0L
  for(sl in sial_letters){
    ## α2-3 ：
    if(any(grepl(paste0("^", sl,"2-[a-z]3$"), list$lin_vec)) ||
       any(grepl(paste0("^[a-z]3-", sl,"2$"), list$lin_vec))){
      cnt23 <- cnt23 + 1L
    }
    ## α2-6 ：
    if(any(grepl(paste0("^", sl,"2-[a-z]6$"), list$lin_vec)) ||
       any(grepl(paste0("^[a-z]6-", sl,"2$"), list$lin_vec))){
      cnt26 <- cnt26 + 1L
    }
  }
  list(cnt23 = cnt23, cnt26 = cnt26)
}

structure_table <- function(list) {
  res_letters <- letters[seq_along(list$seq_raw)]
  seq_raw <- setNames(list$seq_raw, res_letters)
  lin_vec <- list$lin_vec

  edges <- sub("^([a-z])[0-9?]+-([a-z])[0-9?]+$", "\\1-\\2", lin_vec, perl = TRUE)

  descendants <- function(node) {
    kids <- names(father)[father == node]
    if (length(kids) == 0)
      return(character(0))
    c(kids, unlist(lapply(kids, descendants)))
  }
  pairs <- do.call(rbind, strsplit(edges, "-", fixed = TRUE))
  g <- graph_from_edgelist(pairs, directed = FALSE)
  b <- bfs(g, root = "a", father = TRUE)
  father <- setNames(V(g)$name[b$father], V(g)$name)
  father <- father[!is.na(father)]

  branch_man <- names(father)[father == "c"]
  branch_man

  ant_roots <- lapply(branch_man, function(x) names(father)[father == x])
  names(ant_roots) <- branch_man
  # ant_roots are the sugars after branch Mannose
  ant_roots

  ## children of each branch (start from the )
  branch_start <- unlist(ant_roots)
  ant_children <- lapply(branch_start, function(x) descendants(x))
  names(ant_children) <- branch_start
  ant_children

  # Man-β1→4-GlcNAc
  GlcNAc_roots <- seq_raw[branch_start] == "HexNAc"
  ant_cnt <- sum(GlcNAc_roots)
  ant_cnt

  # if is bisecting
  bisecting <- ifelse(any(grepl("^c4-[a-z]1$", lin_vec)), 1, 0)

  # if is complex
  cond1 <- all(GlcNAc_roots)
  # cond2 <- all(sapply(ant_children[GlcNAc_roots], function(x) length(x) >= 1))
  complex <- ifelse(cond1, 1, 0)

  # if is high-mannose
  ant_sugars <- seq_raw[descendants("c")]
  cond1 <- length(ant_sugars) >= 4
  cond2 <- all(ant_sugars %in% c("Mannose", "Galactose"))
  highman <- ifelse(cond1 & cond2, 1, 0)

  # if hybrid
  cond1 <- any(GlcNAc_roots)
  if(length(branch_start) == 0) {
    hybrid <- 0
  } else {
    branch_with_mannose_roots <- seq_raw[branch_start] == "Mannose"
    brancn_with_mannose_desc  <- sapply(branch_start, function(x)
      all(seq_raw[descendants(x)] == "Mannose"))
    cond2 <- any(branch_with_mannose_roots & brancn_with_mannose_desc)
    hybrid <- ifelse(cond1 & cond2, 1, 0)
  }

  # if coreFuc
  corefucosed <- ifelse(any(seq_raw[names(father)[father == "a"]] == "Fucose"), 1, 0)

  list(
    ant_cnt = ant_cnt,
    bisecting = bisecting,
    complex = complex,
    highman = highman,
    hybrid = hybrid,
    corefucosed = corefucosed
  )
}

get_trait <- function(w) {
  ## WURCS --------
  list <- parser_wurcs(w)

  ## Integer
  count_man     <- man_cnt(list)
  count_gal     <- gal_cnt(list)
  count_hexnac  <- hexnac_cnt(list)

  count_sialic  <- sial_cnt(list)
  link_sialic   <- sial_link_cnt(list)
  count_sial_23 <- link_sialic$cnt23
  count_sial_26 <- link_sialic$cnt26
  count_fuc     <- fuc_cnt(list)

  str_table     <- structure_table(list)
  count_antenna <- str_table$ant_cnt
  bisecting     <- str_table$bisecting
  complex       <- str_table$complex
  highman       <- str_table$highman
  hybrid        <- str_table$hybrid
  corefucosed   <- str_table$corefucosed

  res_df <- data.frame(
    Mannose      = count_man,
    Galactose    = count_gal,
    HexNAc       = count_hexnac,
    Sialic       = count_sialic,
    alpha2_3     = count_sial_23,
    alpha2_6     = count_sial_26,
    Fucose       = count_fuc,
    coreFucose   = corefucosed,
    Antennas     = count_antenna,
    Bisecting    = bisecting,
    Complex      = complex,
    HighMan      = highman,
    Hybrid       = hybrid,
    stringsAsFactors = FALSE
  )

  # res_df

  cat(
    "Mannose   :", count_man,     "\n",
    "Galactose :", count_gal,     "\n",
    "HexNAc    :", count_hexnac,  "\n",
    "Sialic    :", count_sialic,  "\n",
    "  ├─α2-3  :", count_sial_23, "\n",
    "  └─α2-6  :", count_sial_26, "\n",
    "Fucose    :", count_fuc,     "\n",
    "  └─core  :", corefucosed,   "\n",
    "Antennas  :", count_antenna, "\n",
    "Bisecting :", bisecting,     "\n",
    "Complex   :", complex,       "\n",
    "HighMan   :", highman,       "\n",
    "Hybrid    :", hybrid,        "\n",
    sep = ""
  )
}

