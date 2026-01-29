#' @keywords internal
"_PACKAGE"

# Selective imports ----

#' @importFrom stats median setNames t.test
#' @importFrom utils read.delim data
#' @importFrom dplyr select group_by summarise bind_rows n all_of
#' @importFrom stringr str_extract_all str_remove
#' @importFrom pbapply pblapply
#' @importFrom car leveneTest
#' @importFrom SummarizedExperiment SummarizedExperiment assays colData rowData
#' @importFrom ggplot2
#'   ggplot aes geom_histogram geom_boxplot geom_jitter
#'   scale_fill_manual labs theme theme_classic
#'   element_text element_blank after_stat
#' @importFrom igraph
#'   graph_from_edgelist graph_from_data_frame neighbors
#'   subcomponent count_subgraph_isomorphisms
#'   bfs V layout_as_tree
#' @importFrom rlang .data
