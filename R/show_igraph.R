library(igraph)

get_igraph_tree <- function(list) {

  # set nodes as a,b,c...
  ids <- c(letters, LETTERS)[seq_along(list$seq_raw)]
  vtab <- data.frame(
    name    = ids,
    residue = list$seq_raw
  )
  names(list$seq_raw) <- ids

  # undirected edge
  undirected_mat <- do.call(rbind, strsplit(list$lin_vec, "-", fixed = TRUE))
  colnames(undirected_mat) <- c("v1", "v2")

  ## create a undirected graph
  g0 <- graph_from_edgelist(undirected_mat, directed = FALSE)

  # from root a to decide direction
  bfs_res <- bfs(g0, root = "a", father = TRUE)
  father  <- setNames(V(g0)$name[bfs_res$father], V(g0)$name)  # child → father
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
  V(g)$is_root <- V(g)$name == "a"


  res_cols <- c(N = "#1f77b4",
                A = "#800080",
                H = "#2ca02c",
                F = "#d62728",
                G = "#7EC1FF")


  V(g)$type <- setNames(list$seq_raw, ids)[V(g)$name]
  V(g)$color <- res_cols[V(g)$type]

  g
}


show_igraph_tree <- function(g) {
  plot(
    g,
    layout         = layout_as_tree(g, root = "a"),
    vertex.label   = V(g)$type,
    vertex.color   = V(g)$color,
    vertex.size    = 15,
    edge.arrow.size = 0.6,
    main           = "Glycan topology coloured by residue type"
  )
}









