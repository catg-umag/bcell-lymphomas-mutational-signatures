suppressMessages(library(factoextra))
library(tidyverse)

# constants
c_palette_mut_pattern <- c("#68b272", "#47b6d7", "#da8a60", "#8d9aca", "#bbc769", "#a45d96")
c_palette_signatures <- c("#7266ae", "#a09ece", "#B4617A", "#48C9B0")
c_palette_signatures_ext <- c("#7266ae", "#A9A7D2", "#c6c4e1", "#bb7087", "#d2a0af", "#5aceb7", "#91decf")
c_grid_style <- element_line(
  size = 0.1, linetype = "solid",
  colour = "#DDDDDD"
)

#' Plots mutational patterns
#'
#' @param data Dataframe containing: $substitution (X>Y), $context (X.Y), and samples (one each column)
#' @param colorby Variable used to color the signature: sample or substitution
plot_patterns_96 <- function(data, colorby = "substitution") {
  colors <- if (colorby == "substitution") {
    c_palette_mut_pattern
  } else {
    c_palette_signatures
  }

  # convert the data to long format
  wl <- tidyr::gather(data, sample, count, -substitution, -context)
  wl$sample <- factor(wl$sample, levels = colnames(data)[-(1:2)]) # match input order

  p <- ggplot(wl) +
    geom_col(aes_string(x = "context", y = "count", fill = colorby), width = 0.75) +
    facet_grid(sample ~ substitution, scales = "free_y") +
    scale_fill_manual(values = colors) +
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

  return(p)
}


plot_signature_contributions <- function(data, from_extraction = TRUE) {
  palette <- if (from_extraction) {
    colorspace::lighten(c_palette_signatures, amount = 0.2)
  } else {
    c_palette_signatures_ext
  }

  p <- ggplot(data = data, aes(x = sample, y = contribution, fill = signature)) +
    geom_col(width = 1) +
    theme_classic() +
    scale_x_discrete(expand = c(0, 1.3)) +
    scale_y_continuous(limits = c(0, 1.001), expand = c(0, 0)) +
    scale_fill_manual(values = palette) +
    theme(
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      text = element_text(size = 9),
      legend.position = "bottom"
    ) +
    ylab("Signature Contribution") +
    xlab("Sample")

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
    xlab("Extracted signatures") +
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
    cex = 0.2,
    k = 3, # Cut in four groups, define in a previous step
    k_colors = c("black", "black", "black"),
    color_labels_by_k = TRUE, # color labels by groups
    rect = TRUE,
    lwd = 0.35,
    rect_border = c("#287091", "#7a4c64", "#CC2936"),
    rect_fill = TRUE,
    show_labels = FALSE,
    labels_track_height = -0.5,
    main = ""
  ) +
    scale_x_continuous(expand = c(0, 1.3))

  return(p)
}