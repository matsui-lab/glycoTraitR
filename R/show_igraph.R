library(igraph)

seq_raw <- list$seq_raw
lin_vec <- list$lin_vec

edge_mat <- do.call(rbind, strsplit(lin_vec, "-"))
g <- graph_from_edgelist(edge_mat, directed = FALSE)

V(g)$letter <- letters[1:length(seq_raw)]
V(g)$type   <- seq_raw

res_cols <- c(N = "#1f77b4",
              F = "#ff7f0e",
              H = "#2ca02c",
              A = "#d62728")
V(g)$color <- res_cols[V(g)$type]


plot(g,
     layout         = layout_as_tree(g, root = "a"),
     vertex.label   = V(g)$type,
     vertex.color   = V(g)$color,
     vertex.size    = 28,
     edge.arrow.size = 0.6,
     main           = "Glycan topology coloured by residue type")
