#' Count monosaccharide residues in a parsed glycan tree
#'
#' Compute residue-level composition traits from a parsed glycan tree.
#' This includes counts of common monosaccharides (Hexose, HexNAc, Neu5Ac,
#' Neu5Gc, Fucose) and the total glycan size.
#'
#' The function operates on glycan trees produced by
#' [`parse_wurcs_structure()`] or [`parse_pGlyco3_structure()`], which
#' represent glycans as residue vectors (`node`) and glycosidic linkages
#' (`edge`).
#'
#' @param tree A parsed glycan tree containing:
#'   * `node`: a vector of residue codes (e.g. `"H"`, `"N"`, `"A"`, `"F"`)
#'   * `edge`: a character vector of parent–child edges
#'
#' @return A named numeric vector containing:
#'   * `GlycanSize` — number of residues
#'   * `Hexose`, `HexNAc`, `Neu5Ac`, `Neu5Gc`, `Fucose`
#'
#' @keywords internal
#' @noRd
count_residues <- function(tree) {
  c(GlycanSize = length(tree$node),
    Hexose = sum(tree$node == "H"),
    HexNAc = sum(tree$node == "N"),
    Neu5Ac = sum(tree$node == "A"),
    Neu5Gc = sum(tree$node == "G"),
    Fucose = sum(tree$node == "F"))
}

#' Construct an igraph representation of a glycan tree
#'
#' Convert a parsed glycan tree (`node` + `edge`) into a directed
#' `igraph` object with parent–child relationships and residue-level
#' metadata suitable for structural motif detection.
#'
#' The resulting graph contains the following vertex attributes:
#' * `name`   — synthetic node label (`"a"`, `"b"`, ...)
#' * `residue` — residue type (H, N, A, F, G)
#' * `type`    — identical to residue (for convenience)
#' * `color`   — color encoding of residue type
#' * `is_root` — TRUE if the vertex is the structural root
#'
#' @param tree A parsed glycan tree with components `node` and `edge`.
#'
#' @return A directed `igraph` object representing the glycan structure.
#'
#' @keywords internal
#' @noRd
build_glycan_igraph <- function(tree) {

  # set nodes as a,b,c...
  ids <- c(letters, LETTERS)[seq_along(tree$node)]
  vtab <- data.frame(
    name    = ids,
    residue = tree$node
  )
  names(tree$node) <- ids

  # undirected edge
  undirected_mat <- do.call(rbind, strsplit(tree$edge, "-", fixed = TRUE))
  colnames(undirected_mat) <- c("v1", "v2")

  ## create a undirected graph
  g0 <- igraph::graph_from_edgelist(undirected_mat, directed = FALSE)

  # from root a to decide direction
  bfs_res <- igraph::bfs(g0, root = "a", father = TRUE)
  father  <- setNames(igraph::V(g0)$name[bfs_res$father], igraph::V(g0)$name)
  father  <- father[!is.na(father)]

  directed_df <- data.frame(
    from = father,                       # parent
    to   = names(father),                # child
    stringsAsFactors = FALSE
  )

  # directed graph
  g <- igraph::graph_from_data_frame(
    d        = directed_df,
    directed = TRUE,
    vertices = vtab
  )
  igraph::V(g)$is_root <- igraph::V(g)$name == "a"

  res_cols <- c("N" = "#1f77b4", "A" = "#800080", "H" = "#2ca02c", "F" = "#d62728", "G" = "#7EC1FF")
  igraph::V(g)$type <- setNames(tree$node, ids)[igraph::V(g)$name]
  igraph::V(g)$color <- res_cols[igraph::V(g)$type]
  g
}

#' Compute structural glycan traits from an igraph glycan tree
#'
#' Evaluate structural glycan features including:
#' * number of antennas
#' * bisecting GlcNAc
#' * complex-type branching
#' * high-mannose characteristics
#' * hybrid-type features
#' * fucosylation (core vs antennary)
#'
#' Additionally, user-defined structural motifs can be quantified
#' using subgraph isomorphism (`igraph::count_subgraph_isomorphisms()`).
#'
#' @details
#' Structural traits are computed from `igraph` objects built by
#' [`build_glycan_igraph()`].
#'
#' User-defined motifs must also be provided as parsed glycan trees
#' (`node` + `edge`), allowing exact structural pattern matching.
#'
#' @param tree A parsed glycan tree.
#' @param motifs Optional named list of user-defined motif trees.
#'
#' @return A named numeric vector combining:
#'   * built-in structural traits (`Antennas`, `Bisect`, `Complex`,
#'     `HighMan`, `Hybrid`, `CoreFuc`, `AntFuc`)
#'   * user-defined motif counts
#'
#' @keywords internal
#' @noRd
compute_structural_traits <- function(tree, motifs) {
  res_letters <- c(letters, LETTERS)[seq_along(tree$node)]
  node <- setNames(tree$node, res_letters)

  # use igraph to represent a glycan structure
  g <- build_glycan_igraph(tree)

  # corefucosed
  kids_of_a <- igraph::neighbors(g, "a", mode = "out")
  corefucosed <- as.integer(any(node[kids_of_a$name] == "F"))

  # antfucosed
  antfucosed <- as.integer((sum(node[kids_of_a$name] == "F") - corefucosed) > 0)

  # find core mannose (N(N("H")))
  coreman <- names(node[node == "H"][1])

  if(!is.na(coreman)){
    # bisecting
    kids_of_coreman <- igraph::neighbors(g, coreman, mode = "out")$name
    cond1 <- length(kids_of_coreman) == 3
    cond2 <- sum(node[kids_of_coreman] == "H") == 2
    cond3 <- sum(node[kids_of_coreman] == "N") == 1
    bisecting <- ifelse(cond1 && cond2 && cond3, 1, 0)

    # number of branches
    branch_man <- node[kids_of_coreman] == "H"
    branch_man <- names(branch_man)[branch_man]
    ant_roots <- lapply(branch_man, function(x)
      igraph::neighbors(g, x, mode = "out")) %>% unlist %>% names()
    ant_cnt <- sum(node[ant_roots] == "N")

    # complex
    cond1 <- ant_cnt >= 2
    cond2 <- all(node[ant_roots] == "N")
    complex <- ifelse(cond1 & cond2, 1, 0)

    # high-mannose
    cond1 <- sum(node == "H") >= 4
    cond2 <- all(node[igraph::subcomponent(g, coreman, mode = "out")] == "H")
    highman <- ifelse(cond1 & cond2, 1, 0)

    # hybrid
    kids_of_each_ant <- lapply(ant_roots, function(x) {
      sub <- igraph::subcomponent(g, x, mode = "out")
      all(node[sub$name] == "H")
    }) %>% unlist
    cond1 <- any(kids_of_each_ant)
    cond2 <- ant_cnt >= 2
    hybrid <- ifelse(cond1 & cond2, 1, 0)
  } else {
    ant_cnt <- 0
    bisecting <- 0
    complex <- 0
    highman <- 0
    hybrid <- 0
  }

  common_traits <- c(
    Antennas = ant_cnt,
    Bisect = bisecting,
    Complex = complex,
    HighMan = highman,
    Hybrid = hybrid,
    CoreFuc = corefucosed,
    AntFuc = antfucosed
  )

  # add user-defined structure
  special_traits <- c()
  if(!is.null(motifs)) {
    n <- length(motifs)
    special_traits <- c()
    for(i in 1:n) {
      trait <- names(motifs)[i]
      g_sub <- build_glycan_igraph(motifs[[i]])

      n_sub <- igraph::count_subgraph_isomorphisms(g_sub, g, method = "vf2",
                                           vertex.color1 = factor(V(g)$type),
                                           vertex.color2 = factor(V(g_sub)$type),
                                           edge.color1 = NULL,
                                           edge.color2 = NULL)
      special_traits[trait] <- n_sub
    }
  }

  c(common_traits, special_traits)
}



