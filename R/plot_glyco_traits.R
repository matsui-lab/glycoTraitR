#' Plot a glycan structure represented as an igraph tree
#'
#' @description
#' Visualize a glycan structure encoded as an igraph object produced by
#' \code{\link{build_glycan_igraph}}. This plot is mainly intended for inspecting
#' glycan topology (branching, residue types, connectivity).
#'
#' @param g
#' An igraph object representing a glycan tree from \code{\link{build_glycan_igraph}}.
#'
#' @return A glycan topology plot is drawn as a side effect.
#'
#'
#' @examples
#' # Example: parse a pGlyco3-style glycan expression into a tree
#' pGlyco_expr <- "(N(N(H(H(H))(H(H)(H)(H(H))))))"
#'
#' # Convert to glycan tree structure
#' tree <- pGlyco3_to_tree(pGlyco_expr)
#'
#' # Explore parsed nodes and edges
#' tree$node
#' tree$edge
#'
#' # Build igraph representation
#' g <- build_glycan_igraph(tree)
#' plot_glycan_tree(g)
#' @export
plot_glycan_tree <- function(g) {
  plot(
    g,
    layout = igraph::layout_as_tree(g, root = "a"),
    vertex.label = igraph::V(g)$type,
    vertex.color = igraph::V(g)$color,
    vertex.size = 30,
    edge.arrow.size = 1,
    main = "Glycan Topology"
  )
}


#' Plot significant glycan trait changes as a heatmap
#'
#' Visualize significant changes in glycan trait abundance and heterogeneity
#' scores across protein- and site-level features.
#'
#' @param res A data frame returned by \code{\link{analyze_hscore_changes}}.
#'
#' @param p_cutoff Numeric p-value cutoff used to select significant
#'   feature-trait pairs. Default is 0.05.
#'
#' @return A \code{ggplot2}/\code{patchwork} object showing significant
#'   feature-trait changes. Tile color indicates the direction and magnitude of
#'   `diff`, and stars indicate significance levels.
#'
#' @details
#' The heatmap displays mean trait abundance (`mu`) and heterogeneity score
#' (`h`) separately. Protein-level results are shown above site-level results
#' when both are available. If either level is absent, only the available panels
#' are shown.
#'
#' @export
plot_hscore_heatmap <- function(res, p_cutoff = 0.05) {

  group_labels <- as.character(res[1, c("group1", "group2")])
  df <- res[res$pval < p_cutoff, ]

  if (nrow(df) == 0) stop("No significant results found.")

  df$sig <- ifelse(df$pval < 0.001, "***", ifelse(df$pval < 0.01, "**", "*"))
  df$plot_diff <- df$diff
  trait_order <- unique(df$trait)

  make_feature_order <- function(level_name) {
    d <- df[df$level == level_name, ]
    if (nrow(d) == 0) return(character(0))

    feature_score <- aggregate(
      abs(diff) ~ feature,
      data = d,
      FUN = max
    )

    feature_score$feature[order(feature_score$`abs(diff)`, decreasing = TRUE)]
  }

  protein_features <- make_feature_order("protein")
  site_features <- make_feature_order("site")

  make_panel <- function(score_type,
                         level_name,
                         features,
                         show_legend = TRUE,
                         show_x_text = TRUE,
                         show_y_axis = TRUE,
                         panel_title = NULL,
                         y_title = NULL) {
    d <- df[
      df$score_type == score_type &
        df$level == level_name,
    ]

    if (nrow(d) == 0 || length(features) == 0 || all(is.na(d$diff))) {
      return(NULL)
    }

    grid <- expand.grid(
      feature = features,
      trait = trait_order,
      stringsAsFactors = FALSE
    )

    d <- merge(
      grid,
      d,
      by = c("feature", "trait"),
      all.x = TRUE
    )

    d$feature <- factor(d$feature, levels = rev(features))
    d$trait <- factor(d$trait, levels = trait_order)

    lim <- max(abs(d$plot_diff), na.rm = TRUE)
    if (!is.finite(lim) || lim == 0) lim <- 1

    ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = trait,
        y = feature,
        fill = plot_diff
      )
    ) +
      ggplot2::geom_tile(
        color = "grey88",
        linewidth = 0.35
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = sig),
        size = 3.3,
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-lim, lim),
        breaks = c(-lim, lim),
        labels = c(
          paste0(group_labels[2], " high"),
          paste0(group_labels[1], " high")
        ),
        na.value = "grey98",
        name = NULL,
        guide = ggplot2::guide_colorbar(
          direction = "horizontal",
          title.position = "top",
          label.position = "bottom",
          ticks = FALSE,
          barwidth = grid::unit(4.5, "cm"),
          barheight = grid::unit(0.35, "cm")
        )
      ) +
      ggplot2::labs(
        title = panel_title,
        x = NULL,
        y = y_title
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position = if (show_legend) "top" else "none",
        legend.justification = "center",

        axis.title.y = if (show_y_axis) {
          ggplot2::element_text(face = "bold", size = 12)
        } else {
          ggplot2::element_blank()
        },

        axis.text.y = if (show_y_axis) {
          ggplot2::element_text(size = 8)
        } else {
          ggplot2::element_blank()
        },

        axis.ticks.y = if (show_y_axis) {
          ggplot2::element_line()
        } else {
          ggplot2::element_blank()
        },

        axis.text.x = if (show_x_text) {
          ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
        } else {
          ggplot2::element_blank()
        },

        axis.ticks.x = if (show_x_text) {
          ggplot2::element_line()
        } else {
          ggplot2::element_blank()
        },

        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        )
      )
  }

  has_protein_mu <- length(protein_features) > 0 &&
    any(df$level == "protein" & df$score_type == "mu")

  has_protein_h <- length(protein_features) > 0 &&
    any(df$level == "protein" & df$score_type == "h")

  has_site_mu <- length(site_features) > 0 &&
    any(df$level == "site" & df$score_type == "mu")

  has_site_h <- length(site_features) > 0 &&
    any(df$level == "site" & df$score_type == "h")

  mu_title_on_protein <- has_protein_mu
  h_title_on_protein  <- has_protein_h


  p_protein_mu <- make_panel(
    score_type = "mu",
    level_name = "protein",
    features = protein_features,
    show_legend = has_protein_mu,
    show_x_text = FALSE,
    show_y_axis = TRUE,
    panel_title = if (mu_title_on_protein) "Mean trait abundance (\u03bc)" else NULL,
    y_title = "Protein level"
  )

  p_protein_h <- make_panel(
    score_type = "h",
    level_name = "protein",
    features = protein_features,
    show_legend = has_protein_h,
    show_x_text = FALSE,
    show_y_axis = FALSE,
    panel_title = if (h_title_on_protein) "Heterogeneity score (h)" else NULL,
    y_title = NULL
  )

  p_site_mu <- make_panel(
    score_type = "mu",
    level_name = "site",
    features = site_features,
    show_legend = !has_protein_mu && has_site_mu,
    show_x_text = TRUE,
    show_y_axis = TRUE,
    panel_title = if (!has_protein_mu && has_site_mu) {
      "Mean trait abundance (\u03bc)"
    } else {
      NULL
    },
    y_title = "Site level"
  )

  p_site_h <- make_panel(
    score_type = "h",
    level_name = "site",
    features = site_features,
    show_legend = !has_protein_h && has_site_h,
    show_x_text = TRUE,
    show_y_axis = FALSE,
    panel_title = if (!has_protein_h && has_site_h) {
      "Heterogeneity score (h)"
    } else {
      NULL
    },
    y_title = NULL
  )

  protein_panels <- Filter(
    Negate(is.null),
    list(p_protein_mu, p_protein_h)
  )

  site_panels <- Filter(
    Negate(is.null),
    list(p_site_mu, p_site_h)
  )

  rows <- list()

  if (length(protein_panels) > 0) {
    rows$protein <- patchwork::wrap_plots(
      protein_panels,
      nrow = 1
    )
  }

  if (length(site_panels) > 0) {
    rows$site <- patchwork::wrap_plots(
      site_panels,
      nrow = 1
    )
  }

  if (length(rows) == 0) {
    stop("No drawable heatmap panels.")
  }

  patchwork::wrap_plots(
    rows,
    ncol = 1,
    heights = if (length(rows) == 2) c(1, 1.2) else 1
  ) +
    patchwork::plot_annotation(
      caption = "* p < 0.05    ** p < 0.01    *** p < 0.001"
    ) &
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        hjust = 0.95,
        vjust = 0.25,
        size = 10
      )
    )
}


#' Plot glycan trait changes as volcano plots
#'
#' Visualize differential glycan trait abundance and heterogeneity results as
#' volcano plots for protein- and site-level features.
#'
#' @param res A data frame returned by \code{\link{analyze_hscore_changes}}
#'
#' @param p_cutoff Numeric p-value cutoff used to define significant points.
#'   Default is 0.05.
#'
#' @param label_size Numeric text size for feature labels. Default is 2.3.
#'
#' @param max_labels Maximum number of significant features to label in each
#'   panel. Features with the smallest p-values are labeled first. Default is 15.
#'
#' @return A \code{ggplot2}/\code{patchwork} object. The x-axis shows group
#'   difference, the y-axis shows `-log10(p-value)`, and colors indicate
#'   significance and direction of change.
#'
#' @export
plot_hscore_volcano <- function(res,
                                p_cutoff = 0.05,
                                label_size = 2.3,
                                max_labels = 15) {

  group_labels <- as.character(res[1, c("group1", "group2")])

  df <- res
  df$xval <- df$diff
  df$neg_log10_p <- -log10(df$pval)

  x_lim <- max(abs(df$xval), na.rm = TRUE)
  if (!is.finite(x_lim) || x_lim == 0) x_lim <- 1

  y_lim <- c(0, max(df$neg_log10_p, na.rm = TRUE) * 1.05)

  df$direction <- "Not significant"
  df$direction[df$pval < p_cutoff & df$xval > 0] <-
    paste0(group_labels[1], " high")
  df$direction[df$pval < p_cutoff & df$xval < 0] <-
    paste0(group_labels[2], " high")

  df$direction <- factor(
    df$direction,
    levels = c(
      "Not significant",
      paste0(group_labels[1], " high"),
      paste0(group_labels[2], " high")
    )
  )

  make_label <- function(trait, feature) {

    peptide <- sub(" \\(.*", "", feature)
    protein <- sub(".*\\(", "", feature)
    protein <- sub("\\)", "", protein)

    paste0(
      trait,
      "\n",
      peptide,
      "\n(",
      protein,
      ")"
    )
  }

  make_panel <- function(score_type,
                         level_name,
                         show_legend = TRUE,
                         show_x_text = TRUE,
                         show_y_axis = TRUE,
                         panel_title = NULL,
                         y_title = NULL) {
    d <- df[
      df$score_type == score_type &
        df$level == level_name,
      ,
      drop = FALSE
    ]

    if (nrow(d) == 0 || all(is.na(d$xval))) {
      return(NULL)
    }

    d_lab <- d[d$pval < p_cutoff, , drop = FALSE]

    if (nrow(d_lab) > max_labels) {
      d_lab <- d_lab[order(d_lab$pval), , drop = FALSE]
      d_lab <- d_lab[seq_len(max_labels), , drop = FALSE]
    }

    d_lab$feature_label <- make_label(
      trait = d_lab$trait,
      feature = d_lab$feature
    )

    ggplot2::ggplot(
      d,
      ggplot2::aes(
        x = xval,
        y = neg_log10_p
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = direction),
        size = 2.1,
        alpha = 0.85
      ) +
      ggplot2::geom_hline(
        yintercept = -log10(p_cutoff),
        linetype = "dashed",
        linewidth = 0.35,
        color = "grey45"
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.35,
        color = "grey45"
      ) +
      ggrepel::geom_text_repel(
        data = d_lab,
        ggplot2::aes(label = feature_label),
        size = label_size,
        max.overlaps = Inf,
        box.padding = 0.6,
        point.padding = 0.35,
        min.segment.length = 0,
        segment.size = 0.25,
        force = 1.8,
        force_pull = 0.15,
        max.time = 2,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(
        values = stats::setNames(
          c("grey75", "#b2182b", "#2166ac"),
          levels(df$direction)
        ),
        drop = FALSE,
        name = NULL
      ) +
      ggplot2::coord_cartesian(
        xlim = c(-x_lim, x_lim),
        ylim = y_lim
      ) +
      ggplot2::labs(
        title = panel_title,
        x = if (show_x_text) "Difference between groups" else NULL,
        y = y_title
      ) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        legend.position = if (show_legend) "top" else "none",
        legend.justification = "center",

        axis.title.y = if (show_y_axis) {
          ggplot2::element_text(face = "bold", size = 12)
        } else {
          ggplot2::element_blank()
        },

        axis.text.y = if (show_y_axis) {
          ggplot2::element_text(size = 9)
        } else {
          ggplot2::element_blank()
        },

        axis.ticks.y = if (show_y_axis) {
          ggplot2::element_line()
        } else {
          ggplot2::element_blank()
        },

        axis.text.x = if (show_x_text) {
          ggplot2::element_text(size = 9)
        } else {
          ggplot2::element_blank()
        },

        axis.ticks.x = if (show_x_text) {
          ggplot2::element_line()
        } else {
          ggplot2::element_blank()
        },

        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5
        )
      )
  }

  has_protein_mu <- any(df$level == "protein" & df$score_type == "mu")
  has_protein_h  <- any(df$level == "protein" & df$score_type == "h")
  has_site_mu    <- any(df$level == "site" & df$score_type == "mu")
  has_site_h     <- any(df$level == "site" & df$score_type == "h")

  p_protein_mu <- make_panel(
    score_type = "mu",
    level_name = "protein",
    show_legend = has_protein_mu,
    show_x_text = TRUE,
    show_y_axis = TRUE,
    panel_title = if (has_protein_mu) "Mean trait abundance (\u03bc)" else NULL,
    y_title = "Protein level"
  )

  p_protein_h <- make_panel(
    score_type = "h",
    level_name = "protein",
    show_legend = has_protein_h,
    show_x_text = TRUE,
    show_y_axis = FALSE,
    panel_title = if (has_protein_h) "Heterogeneity score (h)" else NULL,
    y_title = NULL
  )

  p_site_mu <- make_panel(
    score_type = "mu",
    level_name = "site",
    show_legend = !has_protein_mu && has_site_mu,
    show_x_text = TRUE,
    show_y_axis = TRUE,
    panel_title = if (!has_protein_mu && has_site_mu) {
      "Mean trait abundance (\u03bc)"
    } else {
      NULL
    },
    y_title = "Site level"
  )

  p_site_h <- make_panel(
    score_type = "h",
    level_name = "site",
    show_legend = !has_protein_h && has_site_h,
    show_x_text = TRUE,
    show_y_axis = FALSE,
    panel_title = if (!has_protein_h && has_site_h) {
      "Heterogeneity score (h)"
    } else {
      NULL
    },
    y_title = NULL
  )

  protein_panels <- Filter(
    Negate(is.null),
    list(p_protein_mu, p_protein_h)
  )

  site_panels <- Filter(
    Negate(is.null),
    list(p_site_mu, p_site_h)
  )

  rows <- list()

  if (length(protein_panels) > 0) {
    rows$protein <- patchwork::wrap_plots(
      protein_panels,
      nrow = 1
    )
  }

  if (length(site_panels) > 0) {
    rows$site <- patchwork::wrap_plots(
      site_panels,
      nrow = 1
    )
  }

  if (length(rows) == 0) {
    stop("No drawable volcano panels.")
  }

  patchwork::wrap_plots(
    rows,
    ncol = 1,
    heights = if (length(rows) == 2) c(1, 1.2) else 1
  ) +
    patchwork::plot_annotation(
      caption = paste0(
        "x-axis: group difference; y-axis: -log10(p-value); dashed line: p < ",
        p_cutoff
      )
    ) &
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        hjust = 0.95,
        vjust = 0.25,
        size = 10
      )
    )
}

