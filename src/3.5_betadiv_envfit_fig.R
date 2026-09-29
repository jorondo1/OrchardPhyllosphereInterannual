# PCoA (weighted UniFrac) + envfit taxonomic-class arrows, both barcodes.
# Same ordination and panels as src/2.2_betadiv_fig_main.R, minus the
# ellipses, plus arrows for classes significantly correlated with each ordination.

pacman::p_load(tidyverse, vegan, patchwork, ggrepel, phyloseq, update = FALSE)

source('src/0.0_Config.R')
source('src/utils/beta_div_saver.R') # save_stat_xlsx()

out_dir <- "out/manuscript/supp"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(230726)

betadiv <- readRDS('data/diversity_subsets.rds')
ps <- readRDS('data/ps_objects_full.rds')

# Helper: scores as df, merge metadata (as in 2.2)
score_tibble <- function(pcoa, meta) {
  scores(pcoa, display = "sites") %>%
    as.data.frame() %>%
    rownames_to_column("Sample") %>%
    merge(meta, by = "Sample") %>%
    mutate(
      t_m = factor(
        paste0(Time, ", ", str_to_lower(Management)),
        levels = names(fill_timman)))
}

# Helper: class-level abundance matrix (samples x classes), without
# Unclassified / Incertae sedis classes
class_matrix <- function(sub) {
  agg <- aggregate(t(otu_table(sub)), by = list(class = as.vector(tax_table(sub)[, "Class"])), sum)
  mat <- t(as.matrix(agg[, -1]))
  colnames(mat) <- agg$class
  mat[, !grepl("unclassified|incertae", colnames(mat), ignore.case = TRUE) & colSums(mat) > 0]
}

# Helper: envfit arrows, p <= 0.001 (smallest possible with 999 permutations), length >= 0.4
envfit_arrows <- function(pcoa, sub) {
  cm <- class_matrix(sub)
  stopifnot(identical(rownames(scores(pcoa, display = "sites")), rownames(cm)))
  fit <- envfit(pcoa, cm, permutations = 999)
  as.data.frame(scores(fit, display = "vectors")) %>%
    rownames_to_column("Class") %>%
    mutate(r2 = fit$vectors$r, pval = fit$vectors$pvals, length = sqrt(MDS1^2 + MDS2^2)) %>%
    filter(pval <= 0.001, length >= 0.2)
}

# Bacteria, 2-years -----------------------------------------------------------

pcoa_2y_wuf_16S <- capscale(betadiv$Bacteria$y2$Dist$wuf ~ 1, distance = "wunifrac")

var_2y_wuf_16S <- round(100 * pcoa_2y_wuf_16S$CA$eig / sum(pcoa_2y_wuf_16S$CA$eig), 2)

scores_pcoa_2y_wuf_16S <- score_tibble(pcoa_2y_wuf_16S, betadiv$Bacteria$y2$Meta)

arrows_2y_16S <- envfit_arrows(pcoa_2y_wuf_16S, ps$Bacteria_2y)


# Bacteria, 3-years -----------------------------------------------------------

pcoa_3y_wuf_16S <- capscale(betadiv$Bacteria$y3$Dist$wuf ~ 1, distance = "wunifrac")

var_3y_wuf_16S <- round(100 * pcoa_3y_wuf_16S$CA$eig / sum(pcoa_3y_wuf_16S$CA$eig), 2)

scores_pcoa_3y_wuf_16S <- score_tibble(pcoa_3y_wuf_16S, betadiv$Bacteria$y3$Meta)

arrows_3y_16S <- envfit_arrows(pcoa_3y_wuf_16S, ps$Bacteria_3y)


# Fungi, 2-years ----------------------------------------------------------------

pcoa_2y_wuf_ITS <- capscale(betadiv$Fungi$y2$Dist$wuf ~ 1, distance = "wunifrac")

var_2y_wuf_ITS <- round(100 * pcoa_2y_wuf_ITS$CA$eig / sum(pcoa_2y_wuf_ITS$CA$eig), 2)

scores_pcoa_2y_wuf_ITS <- score_tibble(pcoa_2y_wuf_ITS, betadiv$Fungi$y2$Meta)

arrows_2y_ITS <- envfit_arrows(pcoa_2y_wuf_ITS, ps$Fungi_2y)


# Fungi, 3-years ----------------------------------------------------------------

pcoa_3y_wuf_ITS <- capscale(betadiv$Fungi$y3$Dist$wuf ~ 1, distance = "wunifrac")

var_3y_wuf_ITS <- round(100 * pcoa_3y_wuf_ITS$CA$eig / sum(pcoa_3y_wuf_ITS$CA$eig), 2)

scores_pcoa_3y_wuf_ITS <- score_tibble(pcoa_3y_wuf_ITS, betadiv$Fungi$y3$Meta)

arrows_3y_ITS <- envfit_arrows(pcoa_3y_wuf_ITS, ps$Fungi_3y)

# Plot ---------

# Helper: 2.2's panel (no ellipses) + arrows
pcoa_panel <- function(scores_df, arrows, var_, subtitle) {
  ggplot(scores_df, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
    
    geom_point(shape = 21, stroke = 0, size = 2, alpha = 0.5) +
    geom_segment(
      data = arrows, aes(x = 0, y = 0, xend = MDS1, yend = MDS2), inherit.aes = FALSE,
      arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.3, alpha = 0.8) +
    geom_text_repel(
      data = arrows, inherit.aes = FALSE,
      aes(x = MDS1, y = MDS2, label = Class), 
      force_pull = 1,
      size = 2, segment.color = NA) +
    
    labs(
      x = paste0("PCoA1 ", var_[1], "%"),
      y = paste0("PCoA2 ", var_[2], "%"),
      subtitle = subtitle) +
    
    scale_color_manual(values = color_timman) +
    scale_fill_manual(values = fill_timman)
}

p16S_2y <- pcoa_panel(scores_pcoa_2y_wuf_16S, arrows_2y_16S, var_2y_wuf_16S, "Bacteria - 2-years")
p16S_3y <- pcoa_panel(scores_pcoa_3y_wuf_16S, arrows_3y_16S, var_3y_wuf_16S, "Bacteria - 3-years")
pITS_2y <- pcoa_panel(scores_pcoa_2y_wuf_ITS, arrows_2y_ITS, var_2y_wuf_ITS, "Fungi - 2-years")
pITS_3y <- pcoa_panel(scores_pcoa_3y_wuf_ITS, arrows_3y_ITS, var_3y_wuf_ITS, "Fungi - 3-years")

# Assembled 4-panel figure (Panel A-D) + arrow tables ----------------------------

(p16S_2y + pITS_2y) / (p16S_3y + pITS_3y) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  labs(fill = "Time-Management", color = "Time-Management") &
  theme_pcoa

ggsave(file.path(out_dir, "PCoA_wUF_envfit_panels.pdf"),
       bg = 'white', width = 3000, height = 2300, dpi = 300, units = "px")

save_stat_xlsx(
  file.path(out_dir, "PCoA_wUF_envfit_arrows.xlsx"), 
  rownames_to = NULL, list(
    Bacteria_2y = arrows_2y_16S, Bacteria_3y = arrows_3y_16S,
    Fungi_2y = arrows_2y_ITS, Fungi_3y = arrows_3y_ITS))
