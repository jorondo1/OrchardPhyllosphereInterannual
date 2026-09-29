# Main-figure PCoA panels (weighted UniFrac), both barcodes: Panel A/C
# (Bacteria, 2y/3y), Panel B/D (Fungi, 2y/3y), the combined facet plot, and
# the assembled 4-panel manuscript figure.

pacman::p_load(tidyverse, vegan, patchwork, update = FALSE)

source('src/0.0_Config.R')

dir.create("out/manuscript/betadiv", recursive = TRUE, showWarnings = FALSE)
set.seed(230726)

betadiv <- readRDS('data/diversity_subsets.rds')

# Helper: scores as df, merge metadata
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

# Bacteria, 2-years -----------------------------------------------------------

pcoa_2y_wuf_16S <- capscale(betadiv$Bacteria$y2$Dist$wuf ~ 1, distance = "wunifrac")

var_2y_wuf_16S <- round(100 * pcoa_2y_wuf_16S$CA$eig / sum(pcoa_2y_wuf_16S$CA$eig), 2)

scores_pcoa_2y_wuf_16S <- score_tibble(pcoa_2y_wuf_16S, betadiv$Bacteria$y2$Meta)

p16S_2y <- ggplot(scores_pcoa_2y_wuf_16S, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
  
  geom_point(shape = 21, stroke = 0.2, size = 2) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = .1, linewidth = 0.2) +
  
  labs(
    x = paste0("PCoA1 ", var_2y_wuf_16S[1], "%"), 
    y = paste0("PCoA2 ", var_2y_wuf_16S[2], "%"),
    subtitle = "Bacteria - 2-years") +
  
  scale_color_manual(values = color_timman) +
  scale_fill_manual(values = fill_timman); p16S_2y

# Bacteria, 3-years -----------------------------------------------------------

pcoa_3y_wuf_16S <- capscale(betadiv$Bacteria$y3$Dist$wuf ~ 1, distance = "wunifrac")

var_3y_wuf_16S <- round(100 * pcoa_3y_wuf_16S$CA$eig / sum(pcoa_3y_wuf_16S$CA$eig), 2)

scores_pcoa_3y_wuf_16S <- score_tibble(pcoa_3y_wuf_16S, betadiv$Bacteria$y3$Meta)

p16S_3y <- ggplot(scores_pcoa_3y_wuf_16S, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
  p16S_2y@layers +
  p16S_2y@scales$scales +
  labs(
    x = paste0("PCoA1 ", var_3y_wuf_16S[1], "%"), 
    y = paste0("PCoA2 ", var_3y_wuf_16S[2], "%"),
    subtitle = "Bacteria - 3-years"); p16S_3y

# Fungi, 2-years ----------------------------------------------------------------

pcoa_2y_wuf_ITS <- capscale(betadiv$Fungi$y2$Dist$wuf ~ 1, distance = "wunifrac")

var_2y_wuf_ITS <- round(100 * pcoa_2y_wuf_ITS$CA$eig / sum(pcoa_2y_wuf_ITS$CA$eig), 2)

scores_pcoa_2y_wuf_ITS <- score_tibble(pcoa_2y_wuf_ITS, betadiv$Fungi$y2$Meta)

pITS_2y <- scores_pcoa_2y_wuf_ITS %>% 
  ggplot(aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
  p16S_2y@layers + 
  p16S_2y@scales$scales +
  labs(
    x = paste0("PCoA1 ", var_2y_wuf_ITS[1], "%"),
    y = paste0("PCoA2 ", var_2y_wuf_ITS[2], "%"),
    subtitle = "Fungi - 2-years"); pITS_2y

# Fungi, 3-years ----------------------------------------------------------------

pcoa_3y_wuf_ITS <- capscale(betadiv$Fungi$y3$Dist$wuf ~ 1, distance = "wunifrac")

var_3y_wuf_ITS <- round(100 * pcoa_3y_wuf_ITS$CA$eig / sum(pcoa_3y_wuf_ITS$CA$eig), 2)

scores_pcoa_3y_wuf <- score_tibble(pcoa_3y_wuf_ITS, betadiv$Fungi$y3$Meta)

pITS_3y <- ggplot(scores_pcoa_3y_wuf, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
  p16S_2y@layers + 
  p16S_2y@scales$scales +
  labs(
    x = paste0("PCoA1 ", var_3y_wuf_ITS[1], "%"), 
    y = paste0("PCoA2 ", var_3y_wuf_ITS[2], "%"),
    subtitle = "Fungi - 3-years") ;pITS_3y

# Assembled 4-panel manuscript figure (Panel A-D) -------------------------------

(p16S_2y + pITS_2y) / (p16S_3y + pITS_3y) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  labs(fill = "Time-Management", color = "Time-Management") &
  theme_pcoa

ggsave("out/manuscript/PCoA_wUF_panels.pdf",
       bg = 'white', width = 1700, height = 1300, dpi = 200, units = "px")



# format_var <- function(var_) {
#   var_[rownames(as.data.frame(var_)) %in% c("MDS1", "MDS2")] %>%
#     t() %>% as.data.frame() %>% 
#     rownames_to_column("Variance")
# }
