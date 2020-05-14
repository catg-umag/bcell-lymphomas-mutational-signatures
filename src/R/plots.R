library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
options(scipen = 30)


COLOR_PALETTE <- c("#68b272", "#47b6d7", "#da8a60", "#8d9aca", "#bbc769", "#a45d96")
COLOR2_PALETTE_2 <- c("#E493A5", "#D55672", "#0075A2", "#ADF5FF", "#C5EF77", "#481620", "#C89B7B")
COLOR3_PALETTE_3 <- c("#ADF5FF", "#C5EF77", "#D55672", "#481620", "#0075A2", "#C89B7B")
PALETTE_SPECTRUM <- c("#e5ae7d", COLOR_PALETTE[c(3, 5, 1, 2, 4, 6)])
PALETTE_LYMPHOMA <- c("#CC2936", "#287091")
# PALETTE_LYMPHOMA2 <- c( "#287091", "#287091","#CC2936")
PALETTE_SIGNATURES <- c("#B898FF", "#D4C1FF", "#85C7F2", "#D3D3D3", "#A4D4CA")
PALETTE_MOTIFS <- c("#EF4767", "#1B9AAA", "#DEB841", "#D3D3D3")
PALETTE_FITTING <- c("#294775", "#208eb7", "#99ceeb", "#7a4375", "#946890", "#fd3fbe", "#fd95e8")
paleta_de_choco3 <- c("#d3d3d3", "#3ec0a0", "#3b90bd", "#384abd")
PALETTE_SIGNATURES_NEW <- c(
  "#7266ae", # moradooscuro
  "#a09ece", # moradoclaro
  "#B4617A", # burdeo
  "#48C9B0"
) # verde
PALETTE_LYMPHOMA2 <- c("#4aa5dc", "#4aa5dc", "#1c9262")


fig_save <- function(plot, filename = "plot", formats = c("svg", "pdf", "png"), ...) {
  for (format in formats) {
    ggplot2::ggsave(plot = plot, filename = paste(filename, format, sep = "."), ...)
  }
}

s_theme <- function() {
  t <- theme_classic() + theme(
    legend.position = "none",
    text = element_text(size = 5),
    axis.text.x = element_text(color = "black", size = 7),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7),
    legend.title = element_text(size = 7), legend.text = element_text(size = 5),
    strip.text.x = element_text(size = 7),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 7)
  )

  return(t)
}

s_theme2 <- function() {
  t <- theme_classic() + theme(
    legend.position = "none",
    text = element_text(size = 5),
    axis.text.x = element_text(color = "black", size = 5, angle = 90, hjust = 1),
    axis.text = element_text(size = 5),
    axis.title = element_text(size = 7),
    legend.title = element_text(size = 7), legend.text = element_text(size = 5),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 7)
  )

  return(t)
}


add_significance_annotation <- function(plot, p_value, values) {
  significance <- ifelse(p_value >= 0.01, 0.05,
    ifelse(p_value > 0.0001, 10^ceiling(log10(p_value)), 0.0001)
  )
  max_v <- max(values)
  range_v <- diff(range(values))

  plot <- plot + annotate("text", x = 1.5, y = max_v + range_v * 0.1, label = paste0("p < ", significance), size = 6) +
    annotate("segment", x = 1.2, xend = 1.8, y = max_v + range_v * 0.06, yend = max_v + range_v * 0.06, colour = "black") +
    expand_limits(y = max_v + range_v * 0.12)

  return(plot)
}


plot_mutational_load <- function(data) {
  p <- ggplot(data, aes(x = type, y = n)) +
    geom_boxplot(width = 0.7) +
    stat_boxplot(geom = "errorbar", width = 0.6) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = diff(range(data$n)) / 30, dotsize = 0.7, aes(fill = lymph, color = lymph)) +
    scale_fill_manual(values = PALETTE_LYMPHOMA) +
    scale_color_manual(values = PALETTE_LYMPHOMA) +
    s_theme() +
    labs(y = "Number of Mutations", x = "")

  # statistic test (assumes 2 kinds of lymphomas)
  ttest <- t.test(n ~ lymph, data = data, var.equal = TRUE)
  if (ttest$p.value < 0.05) {
    p <- add_significance_annotation(p, ttest$p.value, data$n)
  }

  return(p)
}

plot_mutational_load_wilcoxon <- function(data) {
  my_comparisons <- list(c("CLL_M", "CLL_UM"), c("CLL_UM", "FL"))
  p <- ggplot(data, aes(x = type, y = n)) +
    geom_boxplot(width = 0.7) +
    stat_boxplot(geom = "errorbar", width = 0.6) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = diff(range(data$n)) / 30, dotsize = 0.7, aes(fill = lymph, color = lymph)) +
    scale_fill_manual(values = PALETTE_LYMPHOMA) +
    scale_color_manual(values = PALETTE_LYMPHOMA) +
    annotate("text", x = 2.25, y = 480, label = paste0("p < ", 0.001), size = 2.25) +
    annotate("segment", x = 1.5, xend = 3, y = 465, yend = 465, colour = "black") +
    annotate("segment", x = 1.5, xend = 1.5, y = 465, yend = 410, colour = "black") +
    annotate("segment", x = 3, xend = 3, y = 465, yend = 455, colour = "black") +
    annotate("segment", x = 1, xend = 2, y = 410, yend = 410, colour = "black") +
    annotate("segment", x = 1, xend = 1, y = 410, yend = 395, colour = "black") +
    annotate("segment", x = 2, xend = 2, y = 410, yend = 395, colour = "black") +
    #          stat_compare_means(comparisons = my_comparisons, label.y=c(440,470), label="p.signif", hide.ns = TRUE, tip.length = c(0.03,0.09),
    #                                      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("p < 0.0001", "p < 0.001", "p < 0.01", "p < 0.05", "ns")))+ # Add pairwise comparisons p-value
    s_theme() +
    labs(y = "Number of Mutations", x = "")

  return(p)
}


plot_vaf <- function(data, vaf_col = "VAF") {
  p <- ggplot(df, aes_string(x = "lymph", y = vaf_col, fill = "lymph")) +
    geom_violin(alpha = 0.5, trim = FALSE, adjust = 0.75, size = 0.1, scale = "count") +
    geom_boxplot(width = .1, outlier.colour = NA, position = position_dodge(width = 0.4)) +
    stat_summary(fun.y = median, geom = "point", size = 6, color = "white") +
    scale_fill_manual(values = PALETTE_LYMPHOMA) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
    s_theme() +
    labs(x = "", y = "Variant Allele Frequency")

  # statistic test (assumes 2 kinds of lymphomas)
  lt <- base::unique(data$lymph)
  mwt <- wilcox.test(df[df$lymph == lt[1], vaf_col], df[df$lymph == lt[2], vaf_col], alternative = "two.sided")
  if (mwt$p.value < 0.05) {
    p <- add_significance_annotation(p, mwt$p.value, data[, vaf_col])
  }

  return(p)
}

plot_vaf_wilcoxon <- function(data) {
  p <- ggplot(df, aes_string(x = "type", y = "tumor_vaf", fill = "type")) +
    geom_violin(alpha = 0.5, trim = FALSE, adjust = 0.9, size = 0.2, scale = "count") +
    geom_boxplot(width = .09, outlier.colour = NA, position = position_dodge(width = 0.4)) +
    stat_summary(fun.y = median, geom = "point", size = 2.5, color = "white") +
    scale_fill_manual(values = PALETTE_LYMPHOMA2) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
    annotate("text", x = 2.25, y = 121, label = paste0("p < ", 0.001), size = 2.25) +
    annotate("segment", x = 1.5, xend = 3, y = 117, yend = 117, colour = "black") +
    annotate("segment", x = 1.5, xend = 1.5, y = 117, yend = 113, colour = "black") +
    annotate("segment", x = 3, xend = 3, y = 117, yend = 108, colour = "black") +
    annotate("segment", x = 1, xend = 2, y = 113, yend = 113, colour = "black") +
    annotate("segment", x = 1, xend = 1, y = 113, yend = 110, colour = "black") +
    annotate("segment", x = 2, xend = 2, y = 113, yend = 110, colour = "black") +
    s_theme() +
    labs(x = "", y = "Variant Allele Frequency")
  return(p)
}


plot_mutation_spectrum <- function(data, percentages = FALSE, CT = FALSE, by = NA, cols = 6) {
  if (CT) {
    colors <- PALETTE_SPECTRUM
    pdata <- data[data$mutation_type != "C>T", ]
  } else {
    colors <- PALETTE_SPECTRUM[-1]
    pdata <- data[!data$mutation_type %in% c("C>T at CpG", "C>T other"), ]
  }
  pdata$sub <- factor(substr(pdata$mutation_type, 1, 3), levels = levels(data$mutation_type))
  if (!is.na(by)) {
    pdata$by <- pdata[, by]
    pdata$nmutations <- sapply(pdata$by, function(x) {
      sum(pdata[pdata$by == x, "n"])
    })
  } else {
    pdata$nmutations <- sum(pdata$n)
  }
  pdata$nmutations <- paste("No. mutations =", as.character(pdata$nmutations))

  vcol <- ifelse(percentages, "perc_m", "n")
  ylabel <- ifelse(percentages, "Relative Contribution (%)", "Mutations")

  p <- ggplot(pdata, aes_string(x = "sub", y = vcol, fill = "mutation_type")) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Substituion Type") +
    labs(x = "", y = ylabel) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.ticks.x = element_blank(),
      strip.background = element_blank()
    )

  if (percentages) {
    pdata$offset <- 0
    if (CT) {
      pdata[pdata$mutation_type == "C>T at CpG", "offset"] <- pdata[pdata$mutation_type == "C>T other", vcol]
    }

    p <- p +
      geom_errorbar(data = pdata, aes(
        x = sub, ymin = perc_m + offset - perc_sd,
        ymax = perc_m + offset + perc_sd
      ), size = 0.3, colour = "black", width = 0.3)
  }

  if (!is.na(by)) {
    p <- p + facet_wrap(by ~ nmutations, ncol = cols)

    nby <- length(unique(pdata$by))
    if (nby > cols && ((nby %% cols != 0) || (nby > 3 * cols))) {
      p <- p + theme(axis.text.x = element_blank(), axis.line.x = element_blank())
    }
  } else {
    p <- p + facet_wrap(~nmutations)
  }

  return(p)
}

plot_mutation_spectrum2 <- function(data, percentages = FALSE, CT = FALSE) {
  if (CT) {
    colors <- PALETTE_SPECTRUM
    pdata <- data[data$mutation_type != "C>T", ]
  } else {
    colors <- PALETTE_SPECTRUM[-1]
    pdata <- data[!data$mutation_type %in% c("C>T at CpG", "C>T other"), ]
  }
  pdata$sub <- factor(substr(pdata$mutation_type, 1, 3), levels = levels(data$mutation_type))

  vcol <- ifelse(percentages, "perc_m", "n")
  ylabel <- ifelse(percentages, "Relative Contribution (%)", "Mutations")

  p <- ggplot(pdata, aes_string(x = "lymph", y = vcol, fill = "mutation_type")) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Substituion Type") +
    labs(x = "", y = ylabel) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.background = element_blank()
    )

  if (percentages) {
    pdata$offset <- 0
    if (CT) {
      pdata[pdata$mutation_type == "C>T at CpG", "offset"] <- pdata[pdata$mutation_type == "C>T other", vcol]
    }

    p <- p +
      geom_errorbar(data = pdata, aes(
        x = lymph, ymin = perc_m + offset - perc_sd,
        ymax = perc_m + offset + perc_sd
      ), size = 0.3, colour = "black", width = 0.3)
  }


  p <- p + facet_wrap(~sub, ncol = 6)

  return(p)
}

plot_patterns <- function(data, colorby = "substitution") {
  colors <- if (colorby == "substitution") {
    COLOR_PALETTE
  } else {
    PALETTE_SIGNATURES_NEW
  }

  # convert the data to long format
  wl <- tidyr::gather(data, sample, count, -substitution, -context)
  wl$sample <- factor(wl$sample, levels = colnames(data)[-(1:2)]) # match input order

  p <- ggplot(wl) +
    geom_col(aes_string(x = "context", y = "count", fill = colorby)) +
    facet_grid(sample ~ substitution, scales = "free_y") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      strip.text = element_text(size = 7),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      legend.position = "none"
    ) +
    xlab("Motif") +
    ylab("Contribution")

  return(p)
}



plot_rainfall <- function(data, genome, chromosomes, by = NA) {
  # dataframe with cummulative sums of chromosomes and centers
  cs_chr <- data.frame(
    chrom = c(chromosomes, NA),
    cumsum = c(0, cumsum(as.numeric(seqlengths(genome)[chromosomes])))
  ) %>%
    dplyr::mutate(diff = c(diff(cumsum), NA)) %>%
    dplyr::mutate(center = cumsum + (diff / 2)) %>%
    dplyr::mutate(chrom = as.character(chrom)) %>%
    dplyr::select(-diff)

  # data for rainfall plot
  data_rf <- data %>%
    dplyr::select(lymph, sample, chrom, pos, substitution) %>%
    dplyr::arrange(lymph, sample, chrom, pos) %>%
    dplyr::mutate(diff = c(-1, diff(pos))) %>%
    dplyr::filter(diff > 0) %>% # discard 1st mutation of each chromosome
    dplyr::mutate(chrom = as.character(chrom)) %>%
    dplyr::left_join(cs_chr, by = "chrom") %>%
    dplyr::mutate(pos_abs = pos + cumsum)

  # plot!
  p <- ggplot(data_rf, aes(x = pos_abs, y = diff)) +
    geom_point(aes(color = substitution), alpha = 0.7, size = 1.5, stroke = 0.3) +
    geom_vline(xintercept = as.vector(cs_chr$cumsum), linetype = "dotted", color = "#888888") +
    scale_color_manual(values = COLOR_PALETTE) +
    theme_bw() +
    scale_y_log10(breaks = c(1e0, 1e2, 1e4, 1e6, 1e8)) +
    scale_x_continuous(label = cs_chr$chr, breaks = cs_chr$center, limits = c(0, max(cs_chr$cumsum)), expand = c(0, 0)) +
    theme(
      text = element_text(size = 12),
      legend.position = "bottom",
      panel.border = element_blank(),
      legend.title = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text.x = element_text(size = 11),
      legend.margin = margin(-20, 0, 0, 0),
      axis.line.x = element_line("#888888", size = 0.3),
      strip.background = element_rect(size = 0, fill = "white")
    ) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(shape = 15, size = 4))) +
    labs(x = "", y = "Genomic Distance")

  if (!is.na(by)) {
    p <- p + facet_wrap(as.formula(paste("~", by)), ncol = 1, scales = "free")
  }

  return(p)
}


plot_aid_motifs <- function(data) {
  colors <- PALETTE_MOTIFS

  p <- ggplot(data, aes("", perc)) +
    geom_col(aes(fill = aid_motif), width = 0.7) +
    scale_fill_manual(values = colors) +
    coord_polar("y") +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.x = element_text(color = "black")
    ) +
    # geom_text(aes(label = paste0(round(perc, 2), "%")), position = position_stack(vjust = 0.5), size = 3) +
    facet_grid(~lymph) +
    labs(title = "Mutations in AID-Motifs", x = "", y = "Percentage [%]", fill = "Motif")

  return(p)
}


plot_mutation_repetition <- function(data) {
  p <- ggplot(data, aes(x = n, y = cs_n, group = lymph)) +
    geom_point(aes(color = lymph)) +
    scale_color_manual(values = PALETTE_LYMPHOMA) +
    geom_path(linetype = "dotted", aes(color = lymph)) +
    theme_classic() +
    theme(text = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 14)) +
    scale_x_continuous(breaks = data$n) +
    labs(y = "Percentage of Total Mutations", x = "Number of Samples", color = "")

  return(p)
}


plot_signature_contribution <- function(data, labels = FALSE, signature_order = c()) {
  # Sort data columns
  data <- data[, sort(colnames(data))]

  # signature order (in bars, without changing names)
  if (length(signature_order) != ncol(data)) {
    signature_order <- colnames(data)
  }
  order_index <- match(signature_order, colnames(data))

  # Normalize by case and set values as percentages
  data <- as.data.frame(data / rowSums(data) * 100)

  # will use this later for sorting the bars
  case_sorted <- data %>%
    tibble::rownames_to_column("case") %>%
    dplyr::mutate(group = as.factor(get_lymph_group(case))) %>%
    dplyr::left_join(
      data.frame(
        group = unique(.$group),
        grmaxsig = colnames(data)[max.col(aggregate(. ~ group, .[, 2:5], sum)[2:4])]
      ),
      by = "group"
    ) %>%
    mutate(grmaxv = apply(., 1, function(x) x[x["grmaxsig"]])) %>%
    # mutate(maxn=max.col(.[2:4]), maxv=apply(.[2:4], 1, max)) %>%
    dplyr::arrange(group, desc(grmaxv)) %>%
    dplyr::mutate(id = 1:nrow(.)) %>%
    dplyr::select(group, case, id)

  # Set groups and transform data to long format
  data <- data %>%
    tibble::rownames_to_column("case") %>%
    tidyr::gather(key = "signature", value = "contribution", -1) %>%
    dplyr::mutate(group = as.factor(get_lymph_group(case))) %>%
    dplyr::select(case, group, everything()) %>%
    as.data.frame()

  # Set a number of 'empty bars' to add at the end of each group
  empty_bar <- 2
  n_signature <- nlevels(as.factor(data$signature))
  to_add <- data.frame(matrix(NA, empty_bar * nlevels(data$group) * n_signature, ncol(data)))
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(data$group), each = empty_bar * n_signature)
  data <- rbind(data, to_add)

  # set ids for ordering the plot
  group_data <- cbind(case_sorted %>% group_by(group) %>% summarize(max = max(id)),
    add = seq(0, (nlevels(data$group) - 1) * empty_bar, empty_bar)
  )

  case_sorted <- case_sorted %>%
    rbind(to_add %>% dplyr::select(group, case) %>% unique() %>% mutate(id = 0)) %>%
    dplyr::left_join(group_data, by = "group") %>%
    dplyr::mutate(id = as.numeric(id), add = as.numeric(add), max = as.numeric(max)) %>%
    dplyr::mutate(id = id + add) %>%
    dplyr::select(-add)
  case_sorted[is.na(case_sorted$case), "id"] <- rowSums(case_sorted[is.na(case_sorted$case), c("id", "max")])

  data <- data %>% dplyr::left_join(case_sorted %>% dplyr::select(-max), by = c("case", "group"))
  for (n in seq(1, empty_bar)) {
    for (gr in unique(data$group)) {
      nmin <- 1 + (n_signature * (n - 1))
      nmax <- n_signature * n
      data[is.na(data$case) & data$group == gr, "id"][nmin:nmax] <-
        data[is.na(data$case) & data$group == gr, "id"][nmin:nmax] + n
    }
  }

  data <- data %>% dplyr::arrange(id)
  data$signature <- factor(data$signature, levels = signature_order)

  # Get the name and the y position of each label
  label_data <- data %>%
    group_by(id, case) %>%
    summarize(tot = sum(contribution))
  nbars <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id - 0.5) / nbars
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)
  label_data$tot <- label_data$tot + 2 # to put som space between the bar and the label

  # prepare a data frame for base lines
  base_data <- data %>%
    group_by(group) %>%
    summarize(start = min(id), end = max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title = mean(c(start, end)))

  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1, ]
  grid_data <- left_join(
    grid_data,
    expand.grid(group = levels(grid_data$group), y = c(0, 25, 50, 75, 100)), # values where the lines are placed
    by = "group"
  )

  # Make the plot
  p <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x = as.factor(id), y = contribution, fill = signature), stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = PALETTE_SIGNATURES_NEW[order_index], breaks = sort(signature_order), name = "Signatures", labels = c("GC", "SBS3+SBS6", "SBS1+SBS5")) +

    # Add 100/75/50/25/0 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = end, y = y, xend = start, yend = y), colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +

    # Add text showing the value of each 100/75/50/25/0 lines
    annotate("text",
      x = rep(max(data$id), 6), y = c(0, 25, 50, 75, 100, 110),
      label = c("0", "25", "50", "75", "100", "%"),
      color = "grey", size = 1.75, angle = 0, fontface = "bold", hjust = .75
    ) +
    ylim(-40, max(label_data$tot, na.rm = T) + 10) + # the +10 y for the "%"
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm"),
      legend.position = c(0.9, 0.23),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 16)
    ) +
    coord_polar(clip = "off") +

    # Add base line information
    geom_segment(
      data = base_data, aes(x = start, y = -5, xend = end, yend = -5),
      colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE
    ) +
    geom_text(
      data = base_data, aes(x = title, y = -10, label = group), hjust = c(0.9, 0.2),
      colour = "black", alpha = 0.8, size = 2, fontface = "bold", inherit.aes = FALSE
    )

  if (labels) {
    p <- p +
      # Add labels on top of each bar
      geom_text(
        data = label_data, aes(x = id, y = tot, label = case, hjust = hjust),
        color = "black", fontface = "bold", alpha = 0.6, size = 1.7, angle = label_data$angle, inherit.aes = FALSE
      ) +
      theme(
        plot.margin = unit(c(-.3, .2, .3, 0.3), "cm"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 6),
        legend.margin = margin(t = -.1, unit = "cm"), legend.key.size = unit(.8, "line")
      ) # plot.margin <- 1st:top,2do:right,3rd:bottom,4to:left
  }

  return(p)
}

plot_signature_contribution_blood <- function(data, labels = FALSE, signature_order = c()) {
  # Sort data columns
  data <- data[, sort(colnames(data))]

  # signature order (in bars, without changing names)
  if (length(signature_order) != ncol(data)) {
    signature_order <- colnames(data)
  }
  order_index <- match(signature_order, colnames(data))

  # Normalize by case and set values as percentages
  data <- as.data.frame(data / rowSums(data) * 100)

  # will use this later for sorting the bars
  case_sorted <- data %>%
    tibble::rownames_to_column("case") %>%
    dplyr::mutate(group = as.factor(get_lymph_group(case))) %>%
    dplyr::left_join(
      data.frame(
        group = unique(.$group),
        grmaxsig = colnames(data)[max.col(aggregate(. ~ group, .[, 2:5], sum)[2:4])]
      ),
      by = "group"
    ) %>%
    mutate(grmaxv = apply(., 1, function(x) x[x["grmaxsig"]])) %>%
    # mutate(maxn=max.col(.[2:4]), maxv=apply(.[2:4], 1, max)) %>%
    dplyr::arrange(group, desc(grmaxv)) %>%
    dplyr::mutate(id = 1:nrow(.)) %>%
    dplyr::select(group, case, id)

  # Set groups and transform data to long format
  data <- data %>%
    tibble::rownames_to_column("case") %>%
    tidyr::gather(key = "signature", value = "contribution", -1) %>%
    dplyr::mutate(group = as.factor(get_lymph_group(case))) %>%
    dplyr::select(case, group, everything()) %>%
    as.data.frame()

  # Set a number of 'empty bars' to add at the end of each group
  empty_bar <- 2
  n_signature <- nlevels(as.factor(data$signature))
  to_add <- data.frame(matrix(NA, empty_bar * nlevels(data$group) * n_signature, ncol(data)))
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(data$group), each = empty_bar * n_signature)
  data <- rbind(data, to_add)

  # set ids for ordering the plot
  group_data <- cbind(case_sorted %>% group_by(group) %>% summarize(max = max(id)),
    add = seq(0, (nlevels(data$group) - 1) * empty_bar, empty_bar)
  )

  case_sorted <- case_sorted %>%
    rbind(to_add %>% dplyr::select(group, case) %>% unique() %>% mutate(id = 0)) %>%
    dplyr::left_join(group_data, by = "group") %>%
    dplyr::mutate(id = as.numeric(id), add = as.numeric(add), max = as.numeric(max)) %>%
    dplyr::mutate(id = id + add) %>%
    dplyr::select(-add)
  case_sorted[is.na(case_sorted$case), "id"] <- rowSums(case_sorted[is.na(case_sorted$case), c("id", "max")])

  data <- data %>% dplyr::left_join(case_sorted %>% dplyr::select(-max), by = c("case", "group"))
  for (n in seq(1, empty_bar)) {
    for (gr in unique(data$group)) {
      nmin <- 1 + (n_signature * (n - 1))
      nmax <- n_signature * n
      data[is.na(data$case) & data$group == gr, "id"][nmin:nmax] <-
        data[is.na(data$case) & data$group == gr, "id"][nmin:nmax] + n
    }
  }

  data <- data %>% dplyr::arrange(id)
  data$signature <- factor(data$signature, levels = signature_order)

  # Get the name and the y position of each label
  label_data <- data %>%
    group_by(id, case) %>%
    summarize(tot = sum(contribution))
  nbars <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id - 0.5) / nbars
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)
  label_data$tot <- label_data$tot + 2 # to put som space between the bar and the label

  # prepare a data frame for base lines
  base_data <- data %>%
    group_by(group) %>%
    summarize(start = min(id), end = max(id) - empty_bar) %>%
    rowwise() %>%
    mutate(title = mean(c(start, end)))

  # prepare a data frame for grid (scales)
  grid_data <- base_data
  grid_data$end <- grid_data$end[c(nrow(grid_data), 1:nrow(grid_data) - 1)] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1, ]
  grid_data <- left_join(
    grid_data,
    expand.grid(group = levels(grid_data$group), y = c(0, 25, 50, 75, 100)), # values where the lines are placed
    by = "group"
  )

  # Make the plot
  p <- ggplot(data) +
    # Add the stacked bar
    geom_bar(aes(x = as.factor(id), y = contribution, fill = signature), stat = "identity", alpha = 0.8) +
    scale_fill_manual(values = PALETTE_SIGNATURES_NEW[order_index], breaks = sort(signature_order), name = "Signatures", labels = c("GC", "SBS3+SBS6", "SBS1+SBS5")) +

    # Add 100/75/50/25/0 lines. I do it at the beginning to make sur barplots are OVER it.
    geom_segment(data = grid_data, aes(x = end, y = y, xend = start, yend = y), colour = "grey", alpha = 1, size = 0.3, inherit.aes = FALSE) +

    # Add text showing the value of each 100/75/50/25/0 lines
    annotate("text",
      x = rep(max(data$id), 6), y = c(0, 25, 50, 75, 100, 110),
      label = c("0", "25", "50", "75", "100", "%"),
      color = "grey", size = 2, angle = 0, fontface = "bold", hjust = .7
    ) +
    ylim(-40, max(label_data$tot, na.rm = T) + 10) + # the +10 y for the "%"
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1, 4), "cm"),
      legend.position = c(0.9, 0.23),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 16)
    ) +
    coord_polar(clip = "off") +

    # Add base line information
    geom_segment(
      data = base_data, aes(x = start, y = -5, xend = end, yend = -5),
      colour = "black", alpha = 0.8, size = 0.6, inherit.aes = FALSE
    ) +
    geom_text(
      data = base_data, aes(x = title, y = -10, label = group), hjust = c(0.9, 0.2),
      colour = "black", alpha = 0.8, size = 2, fontface = "bold", inherit.aes = FALSE
    )

  if (labels) {
    p <- p +
      # Add labels on top of each bar
      geom_text(
        data = label_data, aes(x = id, y = tot, label = case, hjust = hjust),
        color = "black", fontface = "bold", alpha = 0.6, size = 2, angle = label_data$angle, inherit.aes = FALSE
      ) +
      theme(
        plot.margin = unit(c(-.3, .2, .3, 0.5), "cm"), legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 6),
        legend.margin = margin(t = -.01, unit = "cm"), legend.key.size = unit(.8, "line")
      ) # plot.margin <- 1st:top,2do:right,3rd:bottom,4to:left
  }

  return(p)
}

plot_signature_comparison <- function(target_signature, reference_signatures, colorby = "substitution") {
  if (colorby == "substitution") {
    colors <- COLOR_PALETTE
  } else {
    signatures <- colnames(data)[!colnames(data) %in% c("mutation", "substitution", "context")]
    colors <- PALETTE_SIGNATURES_NEW[match(sort(signatures), signatures)]
  }

  wl <- tidyr::gather(reference_signatures, sample, count, -substitution, -context) %>%
    arrange(sample, substitution, context)
  wl$count_target <- rep(target_signature, length(unique(wl$sample)))

  p <- ggplot(wl) +
    geom_col(aes_string(x = "context", y = "count", fill = colorby)) +
    geom_col(aes(x = context, y = count_target), color = "#666666", alpha = 0, size = 0.4) +
    facet_grid(sample ~ substitution) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      strip.text.y = element_text(size = 12),
      strip.text.x = element_text(size = 11),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
      legend.position = "none"
    ) +
    xlab("Motif") +
    ylab("Contribution")

  return(p)
}

plot_fitting <- function(data) {
  p <- ggplot(data = data, aes(x = names, y = contribution, fill = signatures)) +
    geom_col(width = 0.98) + # space between columns
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5), text = element_text(size = 14), legend.position = "bottom", legend.title = element_blank()) +
    labs(x = element_blank(), y = "Contribution") +
    scale_fill_manual(values = PALETTE_FITTING) +
    guides(fill = guide_legend(nrow = 1))

  return(p)
}

plot_mut_by_comp <- function(data, bwidth = 6, size = 0.5) {
  p <- ggplot(data, aes(x = ab, y = n)) +
    geom_boxplot(width = 0.7) +
    facet_grid(~type) +
    geom_signif(comparisons = list(c("A", "B"), c("A", "B")), map_signif_level = T, test = "wilcox.test") +
    stat_boxplot(geom = "errorbar", width = 0.6) +
    geom_dotplot(binaxis = "y", stackdir = "center", binwidth = bwidth, dotsize = size, aes(fill = ab, color = ab)) +
    labs(y = "Number of  AID Mutations by sample")
  return(p)
}

plot_sig_by_comp <- function(data) {
  p <- ggplot(data, aes(x = compartment, y = contrib)) +
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    facet_grid(~type) +
    geom_signif(comparisons = list(c("A", "B")), map_signif_level = T, test = "wilcox.test") +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_dotplot(
      binaxis = "y", stackdir = "center", dotsize = 0.8, alpha = .5, stackratio = 0.5, position = "dodge",
      aes(fill = compartment, color = compartment)
    ) + # binwidth = bwidth, dotsize = size, position = position_dodge(0.001)
    #    scale_y_continuous(limits = c(0, 0.3)) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "white", size = 0)
    ) +
    labs(y = "Sample contribution by compartment")

  return(p)
}

plot_sig_by_type <- function(data) {
  p <- ggplot(data, aes(x = sample, y = contrib)) +
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    stat_boxplot(geom = "errorbar", width = 0.3) +
    geom_dotplot(
      binaxis = "y", stackdir = "center", dotsize = 0.8,
      aes(fill = signature, color = signature)
    ) + # binwidth = bwidth, dotsize = size, position = position_dodge(0.001)
    #    scale_y_continuous(limits = c(0, 0.3)) +
    labs(y = "Signature contribution by disease type")

  return(p)
}

plot_fitting_ab_type <- function(data) {
  p <- ggplot(data, aes(x = compartment, y = contrib, fill = signature)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~type) + # rows = vars(sample)
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(legend.position = "top") +
    #    coord_flip() +
    theme(legend.title = element_blank())
  return(p)
}

plot_fitting_ab <- function(data) {
  p <- ggplot(data, aes(x = sample, y = contrib, fill = signature)) +
    geom_bar(stat = "identity", position = "stack") + # facet_grid(~sample) + #rows = vars(sample)
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(legend.position = "top") +
    #    coord_flip() +
    theme(legend.title = element_blank())

  return(p)
}

plot_fitting_ab_type_perc <- function(data) {
  p <- ggplot(data, aes(x = compartment, y = perc, fill = signatures)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(~type) + # rows = vars(sample)
    s_theme() +
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 10), legend.spacing.x = unit(0.2, "cm"), legend.position = "bottom", legend.box = "vertical") +
    theme(axis.text.x = element_text(size = 12)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(strip.text.x = element_text(size = 12)) +
    theme(
      #    aspect.ratio = 1,
      strip.background = element_blank(),
      strip.placement = "outside"
    ) +
    guides(fill = guide_legend(nrow = 1))
  #    coord_flip() +
  return(p)
}

plot_fitting_ab_sample <- function(data) {
  p <- ggplot(data, aes(x = compartment, y = contrib, fill = signature)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(sample ~ .) + # rows = vars(sample)
    theme(axis.title.y = element_blank()) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text = element_text(size = 10)) +
    theme(axis.text.x = element_text(size = 10)) +
    theme(axis.text.y = element_text(size = 12)) +
    theme(legend.position = "top") +
    coord_flip() +
    theme(legend.title = element_blank())
  return(p)
}

plot_sig_by_comp_dou_fac <- function(data) {
  p <- ggplot(data, aes(x = compartment, y = contribution)) +
    #    geom_boxplot(width = 0.5,outlier.shape= NA, lwd=0.2) +
    scale_color_manual(labels = c("A", "B"), values = PALETTE_LYMPHOMA2[c(2, 3)], name = "Compartment") +
    facet_grid(group ~ signatures, scale = "free") + # type~
    #    geom_signif(comparisons =list(c("A", "B")), map_signif_level = T, test = "t.test", textsize = 3.3,margin_top = 0.05, tip_length=0, vjust=0.2)+ #xmin=x_st, xmax=x_en
    #    compare_means(
    #      contribution ~ compartment,
    #      data,
    #      method = "t.test",
    #      p.adjust.method = "bonferroni") +
    #    stat_compare_means(aes(label = ..p.format..),comparisons = list(c("A", "B")), method = "wilcox.test", tip.length = 0, p.adjust.method = "bonferroni") +
    stat_pvalue_manual(stat.test, label = "p.adj", size = 2) + # vjust =0, vjust =0.18
    #    stat_boxplot(geom = "errorbar", width = 0.2, size=0.3) +
    geom_jitter(aes(fill = compartment, color = compartment), position = position_jitter(0.08), alpha = 1, size = 0.8) +
    stat_summary(
      fun.y = median, fun.ymin = median, fun.ymax = median,
      geom = "crossbar", width = 0.5, size = 0.1
    ) +
    stat_summary(fun.y = median, geom = "point", size = 0.7, color = "black") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 6),
      legend.margin = margin(t = -0.5, unit = "cm")
    ) +
    labs(x = element_blank(), y = "Sample contribution by compartment") +
    guides(fill = FALSE)

  return(p)
}

plot_sig_by_comp_dou_fac2 <- function(data) {
  dat_text <- data.frame(
    label = c("p<0.01", "p<0.01", "p<0.01", "p<0.01"),
    signatures = c("SBS1", "SBS6", "SBS84", "SBS9")
    #    x = c(1.5,1.5,1.5,1.5),
    #    y = c(0.9,0.9,0.9,0.9)
  )

  p <- ggplot(data, aes(x = compartment, y = contribution)) +
    #    geom_boxplot(width = 0.5,outlier.shape= NA, lwd=0.2) +
    scale_color_manual(labels = c("A", "B"), values = PALETTE_LYMPHOMA2[c(2, 3)], name = "Compartment") +
    facet_grid(group ~ signatures, scale = "free") +
    #    annotate("text", x = 'A', y = 0.9 , label = paste0("p < ", 0.01), size = 4) +
    geom_text(data = dat_text, mapping = aes(x = "A", y = 0.8, label = label)) + # , )
    #    annotate("segment", x = 1.5, xend = 3, y = 465, yend = 465, colour = "black") +
    #    annotate("segment", x = 1.5, xend = 1.5, y = 465, yend = 410, colour = "black") +

    #   stat_pvalue_manual(stat.test, label = "p.adj.signif", size = 1.5, vjust =0.18) + #vjust =0
    #    stat_boxplot(geom = "errorbar", width = 0.2, size=0.3) +
    geom_jitter(aes(fill = compartment, color = compartment), position = position_jitter(0.08), alpha = 1, size = 0.8) +
    stat_summary(
      fun.y = median, fun.ymin = median, fun.ymax = median,
      geom = "crossbar", width = 0.5, size = 0.1
    ) +
    stat_summary(fun.y = median, geom = "point", size = 0.7, color = "black") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 6),
      legend.margin = margin(t = -0.5, unit = "cm")
    ) +
    labs(x = element_blank(), y = "Sample contribution by compartment") +
    guides(fill = FALSE)

  return(p)
}

plot_motif_by_sig <- function(data, colores = paleta_de_choco3) {
  p <- ggplot(cosmic_aid_filter, aes_string(x = "Signatures", y = "Percentage", fill = "Motif")) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colores, name = "Motifs") +
    theme_classic() +
    theme(
      text = element_text(size = 7), axis.text.x = element_text(angle = 0, hjust = 0.95),
      legend.position = "bottom",
      #          plot.title = element_text(hjust = 0.5, size = 14),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "white", size = 0),
      legend.key.size = unit(0.5, "line"),
      legend.margin = margin(-5, 0, 0, -12)
    ) +
    coord_flip() +
    guides(fill = guide_legend(reverse = TRUE))

  return(p)
}

plot_dend <- function(data) {
  p <- fviz_dend(res.diana,
    cex = 0.2,
    k = 3, # Cut in four groups, define in a previous step
    k_colors = c("black", "black", "black", "black"),
    color_labels_by_k = TRUE, # color labels by groups
    rect = TRUE,
    lwd = 0.35,
    rect_border = c("#287091", "#7a4c64", "#CC2936"),
    rect_fill = TRUE,
    show_labels = FALSE,
    labels_track_height = -0.5,
    main = ""
  ) +
    scale_x_continuous(expand = expand_scale(add = c(0.3, 0)))

  return(p)
}
