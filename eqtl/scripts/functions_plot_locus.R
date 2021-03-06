# --------------------------------------------------------------
# Make locuscompare or locuszoom plot of GWAS and QTL signals
# - points colored by LD
# Author: Silva Kasela
# --------------------------------------------------------------

library(ggplot2)
library(seqminer)
#source("functions_get_genomatrix.R")

# Simple locuscompare --------

simple_locuscompare <- function(gwas_pval, qtl_pval, gwas_label, qtl_label, main, highlight_index_gwas, highlight_index_qtl, highlight_index_gwas.label = "GWAS hit", highlight_index_qtl.label = "QTL hit") {
  plot(-log10(gwas_pval), -log10(qtl_pval),
     xlab = gwas_label, ylab = qtl_label, main = main,
     ## highlight GWAS and eQTL variant
     col = ifelse(1:length(gwas_pval) == highlight_index_gwas, "orange",
                  ifelse(1:length(gwas_pval) == highlight_index_qtl, "purple", "grey10")),
     pch = ifelse(1:length(gwas_pval) == highlight_index_gwas, 18,
                  ifelse(1:length(gwas_pval) == highlight_index_qtl, 18, 1)),
     cex = ifelse(1:length(gwas_pval) == highlight_index_gwas, 2,
                  ifelse(1:length(gwas_pval) == highlight_index_qtl, 2, 0.65))
     )
  legend.location <- ifelse(-log10(gwas_pval[highlight_index_qtl]) < 2, "topright", "topleft")
  legend(legend.location, c(highlight_index_gwas.label, highlight_index_qtl.label), col = c("orange", "purple"), pch = c(18, 18), bty = "n")
}

# Locuscompare or locuszoom plot with colors ------

locusplot_ld <- function(data, gwas_label, qtl_label, main, pp4,
                         highlight_index_gwas, highlight_index_qtl,
                         highlight_index_gwas.label, highlight_index_qtl.label,
                         highlight_index_gwas.id, highlight_index_qtl.id,
                         ld_variant_id = NULL, ld_variant_group = "qtl",
                         ld_pos_start = NULL, ld_pos_end = NULL,
                         geno = FALSE, genofile = NA,
                         locuscompare = FALSE, locuszoom = FALSE) {
  # dataframe `data` with the at least the following columns: gwas_pval, qtl_pval, phenotype_id
  # matrix `geno`: individuals in rows, variants in columns
  # locuscompare = TRUE: to make locuscompare plot
  # locuscompare = TRUE: to make locuszoom plot
  # highlight_index_xxx.label - variant IDs to show on figures
  # highlight_index_xxx.id - variant IDs to get LD information, needs to match IDs in the geno matrix

  out_fig <- list()

  # Sanity checks - p-values = 0 will be replaced with the most significant p-value that is not 0
  if (sum(data$gwas_pval == 0) > 0) {
    message(unique(data$phenotype_id), " - replacing p-value = 0 with the most significant non-zero p-value")
    find_min <- min(data$gwas_pval[data$gwas_pval != 0])
    data$gwas_pval[data$gwas_pval == 0] <- find_min
  }
  if (sum(data$qtl_pval == 0) > 0) {
    message(unique(data$phenotype_id), " - replacing p-value = 0 with the most significant non-zero p-value")
    find_min <- min(data$qtl_pval[data$qtl_pval != 0])
    data$qtl_pval[data$qtl_pval == 0] <- find_min
  }

  # Figure theme
  theme_set(theme_bw() +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.title = element_text(size = 14),
                  plot.subtitle = element_text(size = 13),
                  axis.title = element_text(size = 12),
                  axis.title.y = element_text(vjust = 2),
                  axis.title.x = element_text(vjust = -0.75),
                  axis.text = element_text(size = 11),
                  strip.text = element_text(size = 11),
                  legend.title = element_text(size = 12),
                  legend.text = element_text(size = 11)))

  # Discrete color gradient from here: https://github.com/tidyverse/ggplot2/issues/2673
  discrete_gradient_pal <- function(colours, bins = 5) {
    ramp <- scales::colour_ramp(colours)
    function(x) {
      if (length(x) == 0) return(character())

      i <- floor(x * bins)
      i <- ifelse(i > bins - 1, bins - 1, i)
      ramp(i/(bins - 1))
    }
  }

  scale_colour_discrete_gradient <- function(..., colours, bins = 5, na.value = "grey50", guide = "colourbar", aesthetics = "colour", colors)  {
    colors <- if (missing(colours))
      colors
    else colours
    continuous_scale(
      aesthetics,
      "discrete_gradient",
      discrete_gradient_pal(colours, bins),
      na.value = na.value,
      guide = guide,
      ...
    )
  }

  # df for highlighting variants on the plot
  ## QTL index variant on the first position, and GWAS index variant on the second position
  highlight_df <- data.frame(qtl_pval = -log10(data$qtl_pval[c(highlight_index_qtl, highlight_index_gwas)]),
                             gwas_pval = -log10(data$gwas_pval[c(highlight_index_qtl, highlight_index_gwas)]),
                             id = c(highlight_index_qtl.id, highlight_index_gwas.id),
                             label = c(highlight_index_qtl.label, highlight_index_gwas.label),
                             pos = data$pos[c(highlight_index_qtl, highlight_index_gwas)],
                             stringsAsFactors = F)

  if (!is.null(ld_variant_id)) {
    # Figure with colors
    if (!is.matrix(geno)) {
      # get genotype matrix
      chr <- unique(data$chr)
      range <- paste0("chr", chr, ":", ld_pos_start, "-", ld_pos_end)
      geno <- get_geno_matrix(genofile = genofile, range = range, variant_id = data$variant_id)
    }
    ld_r <- cor(geno)
    ld_r <- ld_r[ld_variant_id,]
    ld <- ld_r**2
    names(ld) <- sapply(names(ld), function(x){paste(unlist(strsplit(x, "_"))[1:2], collapse = "_")})
    stopifnot(data$id %in% names(ld))
    data$ld <- ld[match(data$id, names(ld))]
    # to change the order of dots plotted
    data <- data[order(data$ld), ]
    rm(ld)

    lz_colors = c("#282973", "#8CCCF0", "#69BD45", "#F9A41A", "#ED1F24")
    # NA color - "#7F7F7F"

    # Highlight the QTL lead SNP as puprle (index SNP for showing LD, if `ld_variant_group` == "qtl") and GWAS lead SNP according to the LD values
    ## QTL index variant on the first position, and GWAS index variant on the second position
    if (ld_variant_group == "qtl") {
      highlight_df$fill <- c("purple",
                             ifelse(data[data$id %in% highlight_index_gwas.id, "ld"] <= 0.2, lz_colors[1],
                                    ifelse(data[data$id %in% highlight_index_gwas.id, "ld"] <= 0.4, lz_colors[2],
                                           ifelse(data[data$id %in% highlight_index_gwas.id, "ld"] <= 0.6, lz_colors[3],
                                                  ifelse(data[data$id %in% highlight_index_gwas.id, "ld"] <= 0.8, lz_colors[4],
                                                         ifelse(data[data$id %in% highlight_index_gwas.id, "ld"] <= 1, lz_colors[5], "grey"))))))
    } else {
      highlight_df$fill <- c(ifelse(data[data$id %in% highlight_index_qtl.id, "ld"] <= 0.2, lz_colors[1],
                                  ifelse(data[data$id %in% highlight_index_qtl.id, "ld"] <= 0.4, lz_colors[2],
                                          ifelse(data[data$id %in% highlight_index_qtl.id, "ld"] <= 0.6, lz_colors[3],
                                                 ifelse(data[data$id %in% highlight_index_qtl.id, "ld"] <= 0.8, lz_colors[4],
                                                        ifelse(data[data$id %in% highlight_index_qtl.id, "ld"] <= 1, lz_colors[5], "grey"))))),
                             "purple")
    }

  }
  if (locuscompare) {
    g <- ggplot(data = data, aes(x = -log10(gwas_pval), y = -log10(qtl_pval))) +
      labs(title = main,
           subtitle = paste0(unique(data$phenotype_id), " locus, PP4 = ", round(100*pp4, 1), "%"),
           x = gwas_label,
           y = qtl_label,
           col = expression(r^2)) +
      scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(.02, .1))) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(.02, .1)))
    if (!is.null(ld_variant_id)) {
      # Figure with colors
      g <- g +
        geom_point(data = data, aes(col = ld), size = 2) +
        scale_colour_discrete_gradient(colours = lz_colors, limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), guide = guide_colourbar(nbin = 100, raster = FALSE, frame.colour = "black", ticks.colour = NA)) +
        geom_label(data = highlight_df, aes(x = gwas_pval, y = qtl_pval, label = label), alpha = 0.25, fill = "grey85", size = 2.5, hjust = "inward", vjust = "bottom", label.padding = unit(1, "mm")) +
        geom_point(data = highlight_df, aes(x = gwas_pval, y = qtl_pval), size = 3, pch = 23, fill = highlight_df$fill, col = "grey25")
      out_fig[["locuscompare"]] <- g
    } else {
    # Figure without colors
    g <- g +
      geom_point(pch = 21, size = 2, fill = "grey75", alpha = 0.7) +
      geom_label(data = highlight_df, aes(x = gwas_pval, y = qtl_pval, label = label), alpha = 0.25, fill = "grey85", size = 2.5, hjust = "inward", vjust = "bottom", label.padding = unit(1, "mm")) +
      geom_point(data = highlight_df, aes(x = gwas_pval, y = qtl_pval), size = 3, pch = 23, fill = "purple")
    out_fig[["locuscompare"]] <- g
    }
  }
  if (locuszoom) {
    # expand scales
 #   scales_df <- data.frame(pval = c(-log10(min(data$gwas_pval)), -log10(min(data$qtl_pval))),
#                            group = c("gwas_pval", "qtl_pval"))
#    scales_df$pval <- ifelse(scales_df$pval < 10, scales_df$pval + 0.75,
#                          ifelse(scales_df$pval < 50, scales_df$pval + 2, scales_df$pval + 5))
#    scales_df$pval <- 10^(-1*scales_df$pval)
    # format data for plotting
    data <- tidyr::pivot_longer(data = data, cols = c("gwas_pval", "qtl_pval"), names_to = "group", values_to = "pval")
    highlight_df <- tidyr::pivot_longer(data = highlight_df, cols = c("gwas_pval", "qtl_pval"), names_to = "group", values_to = "pval")
    group_names <- c(`gwas_pval` = unlist(strsplit(main, " and "))[1],
                     `qtl_pval` = unlist(strsplit(main, " and "))[2])
    # x-axis breaks and labels
    pos_breaks <- pretty(data$pos, n = 4)
    if (nchar(data$pos[1]) > 5) {
      pos <- pos_breaks/1e6
      pos_label <- paste0("Position on chr", data$chr[1], " (Mb)")
    } else {
      pos <- pos_breaks/1e3
      pos_label <- paste0("Position on chr", data$chr[1], " (kb)")
    }

    g <- ggplot(data = data, aes(x = pos, y = -log10(pval))) +
      labs(title = main,
           subtitle = paste0(unique(data$phenotype_id), " locus, PP4 = ", round(100*pp4, 1), "%"),
           x = pos_label,
           y = bquote(-log[10] ~ "(p-value)"),
           col = expression(r^2)) +
      facet_grid(group ~ ., scales = "free_y", labeller = as_labeller(group_names)) +
      scale_x_continuous(breaks = pos_breaks, labels = pos) +
      scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(.02, .1)))

    if (!is.null(ld_variant_id)) {
      # Figure with colors
      g <- g +
        geom_point(aes(col = ld), size = 1.3) +
        scale_colour_discrete_gradient(colours = lz_colors, limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), guide = guide_colourbar(nbin = 100, raster = FALSE, frame.colour = "black", ticks.colour = NA)) +
        geom_label(data = highlight_df, aes(x = pos, y = pval, label = label), alpha = 0.25, fill = "grey85", size = 2.5, hjust = "inward", vjust = "bottom", label.padding = unit(1, "mm")) +
        geom_point(data = highlight_df, aes(x = pos, y = pval), size = 3, pch = 23, fill = highlight_df$fill, col = "grey25")
      out_fig[["locuszoom"]] <- g
    } else {
      # Figure without colors
      g <- g +
        geom_point(pch = 21, size = 1.3, fill = "grey75", alpha = 0.7) +
        geom_label(data = highlight_df, aes(x = pos, y = pval, label = label), alpha = 0.25, fill = "grey85", size = 2.5, hjust = "inward", vjust = "bottom", label.padding = unit(1, "mm")) +
        geom_point(data = highlight_df, aes(x = pos, y = pval), size = 3, pch = 23, fill = "purple")
      out_fig[["locuszoom"]] <- g
    }
  }
  return(out_fig)
}
