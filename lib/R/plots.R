suppressMessages(library(factoextra))
library(tidyverse)

# constants
c_palette_mut_pattern <- c("#68b272", "#47b6d7", "#da8a60", "#8d9aca", "#bbc769", "#a45d96")
c_palette_signatures <- c("#7266ae", "#a09ece", "#b4617a", "#48c9b0", "#cfe7cb")
c_palette_signatures_extended <- c(
  "#7266ae", "#a9a7d2", "#c6c4e1", "#bb7087", "#d2a0af", "#5aceb7", "#91decf", "#cfe7cb", "#e2eaaa"
)
c_palette_aid_motifs <- c("#ef4767", "#1b9aaa", "#deb841", "#d3d3d3")
c_grid_style <- element_line(
  size = 0.1, linetype = "solid",
  colour = "#dddddd"
)

#' Plots mutational patterns
#'
#' @param data Dataframe containing: $substitution (X>Y), $context (X.Y), and samples (one each column)
#' @param colorby Variable used to color the signature: sample or substitution
plot_patterns_96 <- function(data, colorby = "substitution") {
  # convert the data to long format
  wl <- tidyr::gather(data, sample, count, -substitution, -context)
  wl$sample <- factor(wl$sample, levels = colnames(data)[-(1:2)]) # match input order

  p <- ggplot(wl) +
    geom_col(aes_string(x = "context", y = "count", fill = colorby), width = 0.75) +
    facet_grid(sample ~ substitution, scales = "free_y") +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = c_grid_style,
      panel.grid.minor.y = c_grid_style,
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      text = element_text(size = 9),
      legend.position = "none"
    ) +
    xlab("Motif") +
    ylab("Contribution")

  if (colorby == "substitution") {
    p <- p + ggplot2::scale_fill_manual(values = c_palette_mut_pattern)
  } else {
    if (length(unique(wl$sample)) <= length(c_palette_signatures)) {
      p <- p + ggplot2::scale_fill_manual(values = c_palette_signatures)
    } else {
      p <- p + ggplot2::scale_fill_viridis_d()
    }
  }

  return(p)
}


plot_signature_contributions <- function(data, from_extraction = TRUE) {
  nsignatures <- length(unique(data$signature))

  p <- ggplot(data = data, aes(x = sample, y = contribution, fill = signature)) +
    geom_col(width = 1) +
    theme_classic() +
    scale_x_discrete(expand = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1.001), expand = c(0, 0)) +
    theme(
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      text = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    ylab("Signature Contribution") +
    xlab("Sample") +
    guides(fill = guide_legend(nrow = ceiling(nsignatures / 14)))

  if (
    (from_extraction && nsignatures > length(c_palette_signatures)) ||
      (!from_extraction && nsignatures > length(c_palette_signatures_extended))
  ) {
    p <- p + ggplot2::scale_fill_viridis_d()
  } else if (from_extraction) {
    p <- p + ggplot2::scale_fill_manual(values = colorspace::lighten(c_palette_signatures, amount = 0.2))
  } else {
    p <- p + ggplot2::scale_fill_manual(values = c_palette_signatures_extended)
  }

  return(p)
}


plot_extraction_statistics <- function(data) {
  coeff <- ceiling(max(data$mean_cosine_distance) / 0.04) * 0.04
  colors <- c(c_palette_signatures[1], c_palette_signatures[4])

  p <- ggplot(data, aes(x = signatures)) +
    geom_line(aes(y = stability, group = 1), size = 0.75, color = colors[1]) +
    geom_line(aes(y = mean_cosine_distance / coeff, group = 2), size = 0.75, color = colors[2]) +
    scale_y_continuous(
      name = "Signature Average Stability",
      limits = c(0, 1),
      expand = c(0, 0),
      # Add a second axis and specify its features
      sec.axis = sec_axis(~ . * coeff, name = "Mean Sample Cosine Distance", breaks = seq(0, coeff, coeff / 4))
    ) +
    xlab("Extracted Signatures") +
    scale_x_discrete(breaks = data$signatures, expand = c(0.05, 0.05)) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(color = colors[1]),
      axis.title.y.right = element_text(color = colors[2]),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = c_grid_style,
      panel.grid.minor.y = element_blank(),
    )

  return(p)
}


plot_hclust <- function(data) {
  p <- fviz_dend(
    data,
    cex = 0,
    rect = FALSE,
    lwd = 0.35,
    show_labels = FALSE,
    labels_track_height = -0.5,
    main = "",
    ylab = ""
  ) +
    scale_x_continuous(expand = c(0, 1)) +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), plot.margin = unit(c(0, 0, -2, 0), "mm"))

  return(p)
}


plot_aid_motifs <- function(data) {
  p <- ggplot(data, aes(x = name, y = perc)) +
    geom_col(aes(fill = aid_pattern), width = 0.85) +
    theme_classic() +
    coord_flip() +
    scale_fill_manual(values = c_palette_aid_motifs, name = "AID Motif") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    labs(y = "Percentage of Mutations", x = "")

  return(p)
}