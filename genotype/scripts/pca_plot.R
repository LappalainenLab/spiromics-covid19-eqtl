#!/usr/bin/env Rscript
#---------------------------------------------------
# PCA plot using smartpca output files
# Author: Silva Kasela
#---------------------------------------------------

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
#library(RColorBrewer)

# Set theme
theme_set(theme_classic() +
            theme(plot.title = element_text(size = 14),
                  plot.subtitle = element_text(size = 13),
                  axis.title = element_text(size = 13),
                  axis.title.y = element_text(vjust = 2),
                  axis.title.x = element_text(vjust = -0.75),
                  axis.text = element_text(size = 12),
                  strip.text = element_text(size = 12),
                  legend.title = element_text(size = 13),
                  legend.text = element_text(size = 12)))


# Arguments --------------------------------------
args <- commandArgs(trailingOnly = T)
prefix <- args[1]
highlight_samples_file <- args[2]

# Read in data ------------------------------------
evec <- read.table(here("genotype", "pca", "smartpca.pca.evec"), header = F, stringsAsFactors = F)
colnames(evec) <- c("id", paste0("PC", 1:20), "group")
evec$id <- sapply(evec$id, function(x){ unlist(strsplit(x, ":"))[1]})

eval <- read.table(here("genotype", "pca", "smartpca.eval"), header = F, stringsAsFactors = F)
eval$var_explained <- eval$V1/sum(eval$V1)*100
#idx_round <- ifelse(eval$var_explained[5] < 1, 2, 1)
idx_round <- 2

# 1) PCA plot --------------
g <- list()
k1 <- 1
k2 <- 2
for (i in 1:5) {
  evec$x <- evec[, c(k1 + 1)]
  evec$y <- evec[, c(k2 + 1)]
  g[[i]] <- ggplot(data = evec, aes(x = x, y = y)) +
    labs(title = "PCA plot",
         subtitle = prefix,
         x = paste0("PC ", k1, " (", round(eval$var_explained[k1], idx_round), "%)"),
         y = paste0("PC ", k2, " (", round(eval$var_explained[k2], idx_round), "%)")) +
    geom_point(pch = 21)

  k1 <- k1 + 2
  k2 <- k2 + 2
}

# Percent variance explained by top PCs
g[[6]] <- ggplot(data = eval[1:20,], aes(x = 1:20, y = var_explained)) +
  labs(title = paste0("Percent variance explained"),
       subtitle = prefix,
       x = "Genotype PC",
       y = "Variance explained (%)") +
  geom_line() +
  geom_point(pch = 21, fill = "white")
  #scale_x_continuous(breaks = 1:20, labels = 1:20)

# Save the figure
pdf(here("genotype", "pca", "fig", paste0(prefix, ".pca_plot.pdf")), width = 10, height = 7)
g[[1]] + g[[2]] + g[[3]] + g[[4]] + g[[5]] + g[[6]] + plot_layout(ncol = 3)
dev.off()

# 2) PCA plot with highlighted samples ---------------
highlight_samples <- read.table(highlight_samples_file, header = T, sep = "\t", stringsAsFactors = F)
highlight_samples <- highlight_samples[highlight_samples$NWDID_from_master_present_in_freeze9 %in% 1, c("SUBJID", "NWDID_from_master")]

g <- list()
k1 <- 1
k2 <- 2
evec$highlight <- ifelse(evec$id %in% highlight_samples[,2], 1, 0)
evec <- evec[order(evec$highlight),]
for (i in 1:5) {
  evec$x <- evec[, c(k1 + 1)]
  evec$y <- evec[, c(k2 + 1)]
  g[[i]] <- ggplot(data = evec, aes(x = x, y = y, col = as.factor(highlight))) +
    labs(title = "PCA plot",
         subtitle = prefix,
         x = paste0("PC ", k1, " (", round(eval$var_explained[k1], idx_round), "%)"),
         y = paste0("PC ", k2, " (", round(eval$var_explained[k2], idx_round), "%)"),
         col = "w/ RNA-seq") +
    geom_point() + # pch = 21
    scale_color_manual(values = c("grey50", "indianred")) +
    theme(legend.position = "top")

  k1 <- k1 + 2
  k2 <- k2 + 2
}

# Percent variance explained by top PCs
g[[6]] <- ggplot(data = eval[1:20,], aes(x = 1:20, y = var_explained)) +
  labs(title = paste0("Percent variance explained"),
       subtitle = prefix,
       x = "Genotype PC",
       y = "Variance explained (%)") +
  geom_line() +
  geom_point(pch = 21, fill = "white")
  #scale_x_continuous(breaks = 1:20, labels = 1:20)

# Save the figure
pdf(here("genotype", "pca", "fig", paste0(prefix, ".pca_plot.highlight_samples.pdf")), width = 9, height = 7)
g[[1]] + g[[2]] + g[[3]] + g[[4]] + g[[5]] + g[[6]] + plot_layout(ncol = 3)
dev.off()

# Plot of parallel coordinates -----------
evec$x <- NULL
evec$y <- NULL
evec %>%
  pivot_longer(cols = starts_with("PC"), names_to = "variable", values_to = "value") %>%
  mutate("variable" = factor(variable, levels = paste0("PC", 1:20))) %>%
  ggplot(aes(x = variable, y = value, group = id)) +
  labs(title = "Parallel coordinates plot: first 20 PCs",
      subtitle = prefix,
      x = "",
      y = "Value") +
  geom_line(alpha = 0.25) +
  theme(axis.text = element_text(size = 10))
ggsave(here("genotype", "pca", "fig", paste0(prefix, ".parallel_coordinates.pdf")), width = 9, height = 4.5)

## Highlight samples
evec %>%
  filter(id %in% highlight_samples$NWDID_from_master) %>%
  pivot_longer(cols = starts_with("PC"), names_to = "variable", values_to = "value") %>%
  mutate("variable" = factor(variable, levels = paste0("PC", 1:20))) %>%
  ggplot(aes(x = variable, y = value, group = id)) +
  labs(title = "Parallel coordinates plot: first 20 PCs",
       subtitle = paste0(prefix, ": samples with RNA-seq data"),
       x = "",
       y = "Value") +
  geom_line(alpha = 0.25) +
  theme(axis.text = element_text(size = 10))
ggsave(here("genotype", "pca", "fig", paste0(prefix, ".parallel_coordinates.highlight_samples.pdf")), width = 9, height = 4.5)

