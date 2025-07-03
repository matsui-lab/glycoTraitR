parse_glycan <- function(expr) {
  expr <- strsplit(expr, "")[[1]]
  seq_raw <- character()
  lin_vec <- character()
  stack <- c()
  parent <- NULL
  idx <- 1
  for (ch in expr) {
    if (ch == "(") {
      stack <- c(parent, stack)
    } else if (ch == ")") {
      parent <- stack[1]
      stack <- stack[-1]
    } else if (grepl("[A-Za-z]", ch)) {
      label <- letters[idx]
      idx <- idx + 1
      seq_raw <- c(seq_raw, ch)
      if (!is.null(parent)) {
        lin_vec <- c(lin_vec, paste0(parent, "-", label))
      }
      parent <- label
    }
  }
  list(seq_raw = seq_raw, lin_vec = lin_vec)
}



expr <- "(N(N(H(H)(N)(H(H)))))"
list <- parse_glycan(expr); list
expr1 <- build_bracket(list); expr1
expr == expr1


hexose_pGlyco3_cnt  <- function(list) sum(list$seq_raw == "H")
hexnac_pGlyco3_cnt  <- function(list) sum(list$seq_raw == "N")
sialic_pGlyco3_cnt  <- function(list) sum(list$seq_raw == "A")
fucose_pGlyco3_cnt  <- function(list) sum(list$seq_raw == "F")

structure_pGlyco3_table <- function(list) {
  res_letters <- letters[seq_along(list$seq_raw)]
  seq_raw <- setNames(list$seq_raw, res_letters)
  edges <- list$lin_vec

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

  branch_man_children <- lapply(branch_man, function(x) names(father)[father == x])
  names(branch_man_children) <- branch_man
  # ant_roots are the sugars after branch Mannose
  branch_man_children

  ## children of each branch
  branch_start <- unlist(branch_man_children)
  ant_children <- lapply(branch_start, function(x) descendants(x))
  names(ant_children) <- branch_start
  ant_children

  # Hexose-HexNAc
  # branch are defined wherever there is a HexNAc comes next to the Hexose
  ant_roots <- seq_raw[branch_start] == "N"
  ant_cnt <- sum(GlcNAc_roots)
  ant_cnt

  # if is bisecting
  bisecting <- ifelse(length(descendants("c")) == 3, 1, 0)

  # if is complex
  cond1 <- all(ant_roots) & (length(ant_roots) > 1)
  # cond2 <- all(sapply(ant_children[ant_roots], function(x) length(x) >= 1))
  complex <- ifelse(cond1, 1, 0)

  # if is high-mannose
  ant_sugars <- seq_raw[descendants("c")]
  cond1 <- length(ant_sugars) >= 4
  cond2 <- all(ant_sugars %in% c("H"))
  highman <- ifelse(cond1 & cond2, 1, 0)

  # if hybrid
  cond1 <- any(ant_roots)
  if(length(branch_start) == 0) {
    hybrid <- 0
  } else {
    branch_with_mannose_roots <- seq_raw[branch_start] == "H"
    brancn_with_mannose_desc  <- sapply(branch_start, function(x)
      all(seq_raw[descendants(x)] == "H"))
    cond2 <- any(branch_with_mannose_roots & brancn_with_mannose_desc)
    hybrid <- ifelse(cond1 & cond2, 1, 0)
  }

  # if coreFuc
  corefucosed <- ifelse(any(seq_raw[names(father)[father == "a"]] == "F"), 1, 0)

  list(
    ant_cnt = ant_cnt,
    bisecting = bisecting,
    complex = complex,
    highman = highman,
    hybrid = hybrid,
    corefucosed = corefucosed
  )
}

get_trait <- function(expr) {
  ## pGlyco3 --------
  list <- parse_glycan(expr)

  ## Integer
  count_hexose  <- hexose_pGlyco3_cnt(list)
  count_hexnac  <- hexnac_pGlyco3_cnt(list)
  count_sialic  <- sialic_pGlyco3_cnt(list)
  count_fucose  <- fucose_pGlyco3_cnt(list)
  count_hexose;
  count_hexnac;
  count_sialic;
  count_fucose;

  str_table     <- structure_pGlyco3_table(list)
  bisecting     <- str_table$bisecting
  count_antenna <- str_table$ant_cnt
  complex       <- str_table$complex
  highman       <- str_table$highman
  hybrid        <- str_table$hybrid
  corefucosed   <- str_table$corefucosed

  # res_df <- data.frame(
  #   Hexose       = count_hexose,
  #   HexNAc       = count_hexnac,
  #   Sialic       = count_sialic,
  #   Fucose       = count_fucose,
  #   coreFucose   = corefucosed,
  #   Antennas     = count_antenna,
  #   Bisecting    = bisecting,
  #   Complex      = complex,
  #   HighMan      = highman,
  #   Hybrid       = hybrid,
  #   stringsAsFactors = FALSE
  # )

  # res_df

  cat(
    "Hexose    :", count_hexose,  "\n",
    "HexNAc    :", count_hexnac,  "\n",
    "Sialic    :", count_sialic,  "\n",
    "Fucose    :", count_fucose,  "\n",
    "  └─core  :", corefucosed,   "\n",
    "Antennas  :", count_antenna, "\n",
    "Bisecting :", bisecting,     "\n",
    "Complex   :", complex,       "\n",
    "HighMan   :", highman,       "\n",
    "Hybrid    :", hybrid,        "\n",
    sep = ""
  )
}
