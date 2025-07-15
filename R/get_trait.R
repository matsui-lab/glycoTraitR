hexose_cnt  <- function(list) sum(list$seq_raw == "H")
hexnac_cnt  <- function(list) sum(list$seq_raw == "N")
neu5ac_cnt  <- function(list) sum(list$seq_raw == "A")
neu5gc_cnt  <- function(list) sum(list$seq_raw == "G")
fucose_cnt  <- function(list) sum(list$seq_raw == "F")

structure_table <- function(list) {
  res_letters <- c(letters, LETTERS)[seq_along(list$seq_raw)]
  seq_raw <- setNames(list$seq_raw, res_letters)
  edges <- list$lin_vec

  # use igraph to represent a glycan structure
  pairs <- do.call(rbind, strsplit(edges, "-"))
  g <- graph_from_edgelist(pairs, directed = TRUE)
  # subcomponent(g, v = "c", mode = "out")

  # corefucosed
  kids_of_a <- neighbors(g, "a", mode = "out")
  corefucosed <- as.integer(any(seq_raw[kids_of_a$name] == "F"))

  # find core mannose
  coreman <- names(seq_raw[seq_raw == "H"][1])

  # bisecting
  kids_of_coreman <- neighbors(g, coreman, mode = "out")$name
  cond1 <- length(kids_of_coreman) == 3
  cond2 <- sum(seq_raw[kids_of_coreman] == "H") == 2
  cond3 <- sum(seq_raw[kids_of_coreman] == "N") == 1
  bisecting <- ifelse(cond1 && cond2 && cond3, 1, 0)

  # number of branches
  branch_man <- seq_raw[kids_of_coreman] == "H"
  branch_man <- names(branch_man)[branch_man]
  ant_roots <- lapply(branch_man, function(x) neighbors(g, x, mode = "out")) %>% unlist %>% names()
  ant_cnt <- sum(seq_raw[ant_roots] == "N")

  # complex
  cond1 <- ant_cnt >= 2
  cond2 <- all(seq_raw[ant_roots] == "N")
  complex <- ifelse(cond1 & cond2, 1, 0)

  # high-mannose
  cond1 <- sum(seq_raw == "H") >= 4
  cond2 <- all(seq_raw[subcomponent(g, coreman, mode = "out")] == "H")
  highman <- ifelse(cond1 & cond2, 1, 0)

  # hybrid
  kids_of_each_ant <- lapply(ant_roots, function(x) {
    sub <- subcomponent(g, x, mode = "out")
    # sub <- as.numeric(sub)
    all(seq_raw[sub$name] == "H")
    }) %>% unlist
  cond1 <- any(kids_of_each_ant)
  cond2 <- ant_cnt >= 2
  hybrid <- ifelse(cond1 & cond2, 1, 0)

  list(
    ant_cnt = ant_cnt,
    bisecting = bisecting,
    complex = complex,
    highman = highman,
    hybrid = hybrid,
    corefucosed = corefucosed
  )
}

get_trait <- function(list) {

  count_N  <- hexose_cnt(list)
  count_H  <- hexnac_cnt(list)
  count_A  <- neu5ac_cnt(list)
  count_G  <- neu5gc_cnt(list)
  count_F  <- fucose_cnt(list)

  topo      <- structure_table(list)
  ant_cnt   <- topo$ant_cnt
  bisecting <- topo$bisecting
  complex   <- topo$complex
  highman   <- topo$highman
  hybrid    <- topo$hybrid
  core_fuc  <- topo$corefucosed

  # cat(
  #   "Hexose    : ", count_N,  "\n",
  #   "HexNAc    : ", count_H,  "\n",
  #   "Neu5Ac    : ", count_A,  "\n",
  #   "Neu5Gc    : ", count_G,  "\n",
  #   "Fucose    : ", count_F,  "\n",
  #   "  └─core  : ", core_fuc, "\n",
  #   "Antennas  : ", ant_cnt,  "\n",
  #   "Bisecting : ", bisecting, "\n",
  #   "Complex   : ", complex,   "\n",
  #   "HighMan   : ", highman,   "\n",
  #   "Hybrid    : ", hybrid,    "\n",
  #   sep = "")

  list(
    Hexose   = count_N,
    HexNAc   = count_H,
    Neu5Ac   = count_A,
    Neu5Gc   = count_G,
    Fucose   = count_F,
    coreFuc  = core_fuc,
    Antennas = ant_cnt,
    Bisect   = bisecting,
    Complex  = complex,
    HighMan  = highman,
    Hybrid   = hybrid
  )
}

get_trait_from_wurcs <- function(w) {
  w %>% wurcs_to_tree() %>% get_trait() %>% unlist()
}

