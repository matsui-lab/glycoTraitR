cnt_table <- function(list) {
  c(GlycanSize = length(list$node),
    Hexose = sum(list$node == "H"),
    HexNAc = sum(list$node == "N"),
    Neu5Ac = sum(list$node == "A"),
    Neu5Gc = sum(list$node == "G"),
    Fucose = sum(list$node == "F"))
}

get_igraph_tree <- function(list) {

  # set nodes as a,b,c...
  ids <- c(letters, LETTERS)[seq_along(list$node)]
  vtab <- data.frame(
    name    = ids,
    residue = list$node
  )
  names(list$node) <- ids

  # undirected edge
  undirected_mat <- do.call(rbind, strsplit(list$edge, "-", fixed = TRUE))
  colnames(undirected_mat) <- c("v1", "v2")

  ## create a undirected graph
  g0 <- igraph::graph_from_edgelist(undirected_mat, directed = FALSE)

  # from root a to decide direction
  bfs_res <- igraph::bfs(g0, root = "a", father = TRUE)
  father  <- setNames(igraph::V(g0)$name[bfs_res$father], igraph::V(g0)$name)  # child → father
  father  <- father[!is.na(father)]

  directed_df <- data.frame(
    from = father,                       # parent
    to   = names(father),                # child
    stringsAsFactors = FALSE
  )

  # directed graph
  g <- graph_from_data_frame(
    d        = directed_df,
    directed = TRUE,
    vertices = vtab
  )
  igraph::V(g)$is_root <- igraph::V(g)$name == "a"

  res_cols <- c("N" = "#1f77b4", "A" = "#800080", "H" = "#2ca02c", "F" = "#d62728", "G" = "#7EC1FF")
  igraph::V(g)$type <- setNames(list$node, ids)[igraph::V(g)$name]
  igraph::V(g)$color <- res_cols[igraph::V(g)$type]
  g
}


struct_table <- function(list, user_defined_traits) {
  res_letters <- c(letters, LETTERS)[seq_along(list$node)]
  node <- setNames(list$node, res_letters)
  edges <- list$edge

  # use igraph to represent a glycan structure
  g <- get_igraph_tree(list)

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
  if(!is.null(user_defined_traits)) {
    n <- length(user_defined_traits)
    special_traits <- c()
    for(i in 1:n) {
      trait <- names(user_defined_traits)[i]
      g_sub <- get_igraph_tree(user_defined_traits[[i]])

      n_sub <- count_subgraph_isomorphisms(g_sub, g, method = "vf2",
                                           vertex.color1 = factor(V(g)$type),
                                           vertex.color2 = factor(V(g_sub)$type),
                                           edge.color1 = NULL,
                                           edge.color2 = NULL)
      special_traits[trait] <- n_sub
    }
  }

  c(common_traits, special_traits)
}



