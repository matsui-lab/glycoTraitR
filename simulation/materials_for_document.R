library(pbapply)
decipher_input <- readRDS("./simulation/decipher_input.rds")
pGlyco3_input <- readRDS("./simulation/pGlyco3_input.rds")
meta <- readRDS("./simulation/meta.rds")

w <- decipher_input$GlycanStructure
ws <- lapply(w, wurcs_to_tree)

st <- pGlyco3_input$GlycanStructure
sts <- lapply(st, pGlyco3_to_tree)

# Show trait definition
sts[[1]] %>% get_igraph_tree() %>% show_igraph_tree()
ws[[1]] %>% get_igraph_tree() %>% show_igraph_tree()

#### Example 1 #####
list <- sts[[633]]
g <- get_igraph_tree(list)

list_sub <- list()
list_sub$node <- c("H", "H", "H")
list_sub$edge <- c("a-b", "b-c")
g_sub <- get_igraph_tree(list_sub)

n_sub <- count_subgraph_isomorphisms(g_sub, g, method = "vf2",
                                     vertex.color1 = factor(V(g)$type),
                                     vertex.color2 = factor(V(g_sub)$type),
                                     edge.color1 = NULL,
                                     edge.color2 = NULL)
par(mfrow=c(1,2))
show_igraph_tree(g_sub); show_igraph_tree(g); n_sub
#### Example 2 #####
list <- sts[[2478]]
g <- get_igraph_tree(list)

list_sub <- list()
list_sub$node <- c("H", "N", "F")
list_sub$edge <- c("a-b", "b-c")
g_sub <- get_igraph_tree(list_sub)

n_sub <- count_subgraph_isomorphisms(g_sub, g, method = "vf2",
                                     vertex.color1 = factor(V(g)$type),
                                     vertex.color2 = factor(V(g_sub)$type),
                                     edge.color1 = NULL,
                                     edge.color2 = NULL)
par(mfrow=c(1,2))
show_igraph_tree(g_sub); show_igraph_tree(g); n_sub

################################################################################
# User can define their only traits
user_defined_traits <- list(
  "Sub_struct_1" = list(node = c("H", "N", "F"), edge = c("a-b", "b-c")),
  "Sub_struct_2" = list(node = c("H", "H", "H"), edge = c("a-b", "b-c"))
)
list <- ws[[2478]]
get_trait(list, user_defined_traits)

