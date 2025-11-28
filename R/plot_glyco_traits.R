show_igraph_tree <- function(g) {
  plot(
    g,
    layout         = igraph::layout_as_tree(g, root = "a"),
    vertex.label   = igraph::V(g)$type,
    vertex.color   = igraph::V(g)$color,
    vertex.size    = 30,
    edge.arrow.size = 1,
    main           = "Glycan Topology"
  )
}

get_plots <- function(trait_se, group_list, trait, level) {
  assays <- assays(trait_se)
  coldata <- colData(trait_se)
  rowdata <- rowData(trait_se)

  abund <- assays[[trait]][level, ]
  label <- coldata[[names(group_list)]]

  # select psm corresponding to files in the group
  group     <- group_list[[names(group_list)]]
  col_ind   <- label %in% group
  label_new <- label[col_ind]
  abund_new <- abund[col_ind]

  # select non NA psm
  ind <- !is.na(abund_new)
  label_new <- label_new[ind]
  abund_new <- abund_new[ind]

  df_plot <- data.frame(Abundance = as.numeric(abund_new),
                        Group     = factor(label_new, levels = group))

  ## freq hist
  p_freq <-
    ggplot(df_plot, aes(x = Abundance, fill = Group)) +
    geom_histogram(
      aes(y = ..count..),
      position = "dodge",
      color = "black",
      binwidth = max(1, ceiling(max(df_plot$Abundance) / 20))
    ) +
    scale_fill_manual(values = setNames(c("#e31a1c", "#1f78b4"), group)) +
    labs(x = paste(trait, "Count", "at", level), y = "PSM count") +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14)
    )

  p_box <- ggplot(df_plot, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(width = 0.5, show.legend = FALSE) +
    geom_jitter(size = 1, width = 0.1, show.legend = FALSE) +
    scale_fill_manual(values = setNames(c("#e31a1c", "#1f78b4"), group)) +
    labs(x = NULL, y = paste(trait, "Count", "at", level)) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14)
    )

  return(list(p_hist = p_freq, p_box = p_box))
}


