#' Construct an igraph representation of a glycan tree
#'
#' Convert a parsed glycan tree into a directed
#' `igraph` object with parent–child relationships and residue-level
#' metadata suitable for structural motif detection.
#'
#' The resulting graph contains the following vertex attributes:
#' \itemize{
#'   \item \code{name}    — synthetic node label (`"a"`, `"b"`, ...)
#'   \item \code{residue} — residue type (H, N, A, F, G)
#'   \item \code{type}    — identical to residue (for convenience)
#'   \item \code{color}   — color encoding of residue type
#'   \item \code{is_root} — TRUE if the vertex is the structural root
#' }
#'
#' @param tree A parsed glycan tree from \code{\link{pGlyco3_to_tree}} or \code{\link{wurcs_to_tree}}.
#'
#' @return A directed `igraph` object representing the glycan structure.
#'
#' @examples
#' # Example: parse a pGlyco3 monosaccharide expression into a glycan tree
#' pGlyco_expr <- "(N(N(H(H(H))(H(H)(H)(H(H))))))"
#'
#' # Convert expression into a parsed tree structure
#' tree <- pGlyco3_to_tree(pGlyco_expr)
#' g <- build_glycan_igraph(tree)
#' g
#'
#' @export
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
  father <- setNames(igraph::V(g0)$name[bfs_res$father], igraph::V(g0)$name)
  father <- father[!is.na(father)]

  directed_df <- data.frame(
    from = father,
    to = names(father),
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

#' Count monosaccharide residues in a parsed glycan tree
#'
#' Compute residue-level composition traits from a parsed glycan tree.
#' This includes counts of common monosaccharides (Hexose, HexNAc, Neu5Ac,
#' Neu5Gc, Fucose) and the total glycan size.
#'
#' @details
#' The function operates on glycan trees produced by
#' \code{\link{build_glycan_tree}} or \code{\link{parse_pGlyco3_structure}}, which
#' represent glycans as residue vectors (`node`) and glycosidic linkages
#' (`edge`).
#'
#' @param tree A parsed glycan tree from \code{\link{pGlyco3_to_tree}} or \code{\link{wurcs_to_tree}}.
#'
#' @return A named numeric vector containing:
#'   `GlycanSize`, `Hexose`, `HexNAc`, `Neu5Ac`, `Neu5Gc`, `Fucose`
#'
#' @keywords internal
#' @noRd
count_residues <- function(tree) {
  c(
    GlycanSize = length(tree$node),
    Hexose = sum(tree$node == "H"),
    HexNAc = sum(tree$node == "N"),
    Neu5Ac = sum(tree$node == "A"),
    Neu5Gc = sum(tree$node == "G"),
    Fucose = sum(tree$node == "F")
  )
}

#' Compute structural glycan traits from an igraph glycan tree
#'
#' Evaluate build-in structural glycan traits including:
#' \itemize{
#'   \item \code{Antenna numbers}
#'   \item \code{Bisecting-type}
#'   \item \code{Complex-type}
#'   \item \code{High-mannose-type}
#'   \item \code{Hybrid-type}
#'   \item \code{Core-fucosed}
#'   \item \code{Antennary-fucosed}
#' }
#'
#' @param tree A parsed glycan tree from \code{\link{pGlyco3_to_tree}} or \code{\link{wurcs_to_tree}}.
#'
#' @return A named numeric vector of built-in structural traits and
#' user-defined motif counts
#'
#' @keywords internal
#' @noRd
compute_structural_traits <- function(tree) {
  res_letters <- c(letters, LETTERS)[seq_along(tree$node)]
  node <- setNames(tree$node, res_letters)
  # use igraph to represent a glycan structure
  g <- build_glycan_igraph(tree)
  # corefucosed
  kids_of_a <- igraph::neighbors(g, "a", mode = "out")
  corefucosed <- as.integer(any(node[kids_of_a$name] == "F"))
  # antfucosed
  antfucosed <- as.integer((sum(node[kids_of_a$name] == "F") - corefucosed) > 0)

  coreman <- names(node[node == "H"][1])
  if (!is.na(coreman)) {
    # bisecting
    kids_of_coreman <- igraph::neighbors(g, coreman, mode = "out")$name
    cond1 <- length(kids_of_coreman) == 3
    cond2 <- sum(node[kids_of_coreman] == "H") == 2
    cond3 <- sum(node[kids_of_coreman] == "N") == 1
    bisecting <- ifelse(cond1 && cond2 && cond3, 1, 0)

    # number of branches
    branch_man <- node[kids_of_coreman] == "H"
    branch_man <- names(branch_man)[branch_man]
    ant_roots <- lapply(branch_man, function(x) {
      igraph::neighbors(g, x, mode = "out")
    }) %>%
      unlist() %>%
      names()
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
    }) %>% unlist()
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
  c(
    Antennas = ant_cnt, Bisect = bisecting, Complex = complex,
    HighMan = highman, Hybrid = hybrid, CoreFuc = corefucosed, AntFuc = antfucosed
  )
}


#' Compute user-defined structural glycan traits via subgraph isomorphism
#'
#' Evaluate custom glycan structural motifs within a parsed glycan tree.
#' Each user-defined motif is represented as a parsed tree (`node` + `edge`),
#' converted into an `igraph` structure, and matched against the full glycan
#' using subgraph isomorphism (`igraph::count_subgraph_isomorphisms()`).
#'
#' @details
#' This function allows end users to define their own glycan structural motifs—
#' such as specific antenna patterns, branch-residue
#' combinations—and quantify how often these motifs appear in a glycan.
#' The output is a named numeric vector where names correspond to motif names.
#'
#' @param tree A parsed glycan tree from \code{\link{pGlyco3_to_tree}} or \code{\link{wurcs_to_tree}}.
#'
#' @param motifs Optional named list of user-defined motif trees.
#'   Each motif must include:
#'   \itemize{
#'     \item `node`: character vector of residue codes
#'     \item `edge`: edges describing motif topology
#'   }
#'
#' @return
#' A named numeric vector giving the count of each user-defined motif.
#' Returns an empty vector if \code{motifs = NULL}.
#'
#' @keywords internal
#' @noRd
compute_userdefined_traits <- function(tree, motifs) {
  res_letters <- c(letters, LETTERS)[seq_along(tree$node)]
  node <- setNames(tree$node, res_letters)

  # use igraph to represent a glycan structure
  g <- build_glycan_igraph(tree)

  # add user-defined structure
  ud_traits <- c()
  if (!is.null(motifs)) {
    n <- length(motifs)
    special_traits <- c()
    for (i in seq_len(n)) {
      trait <- names(motifs)[i]
      g_sub <- build_glycan_igraph(motifs[[i]])
      n_sub <- igraph::count_subgraph_isomorphisms(g_sub, g,
        method = "vf2",
        vertex.color1 = factor(V(g)$type),
        vertex.color2 = factor(V(g_sub)$type),
        edge.color1 = NULL,
        edge.color2 = NULL
      )
      ud_traits[trait] <- n_sub
    }
  }

  ud_traits
}
