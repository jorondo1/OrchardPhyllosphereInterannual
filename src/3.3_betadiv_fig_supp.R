# Supp companion to 2.2's main wUF figure: essentially the same 4-panel
# PCoA figure (Panel A/C Bacteria 2y/3y, Panel B/D Fungi 2y/3y, combined
# facet plot, assembled patchwork), but on Bray-Curtis distances.

pacman::p_load(tidyverse, vegan, patchwork, update = FALSE)

source('src/0.0_Config.R')
betadiv <- readRDS('data/diversity_subsets.rds')

dir.create("out/manuscript/betadiv", recursive = TRUE, showWarnings = FALSE)
set.seed(230726)

# Bacteria, 2-years -----------------------------------------------------------

pcoa_2y_bray <- capscale(betadiv$Bacteria$y2$Dist$bray ~ 1, distance = "bray")

variance_2y_bray <- round(100 * pcoa_2y_bray$CA$eig / sum(pcoa_2y_bray$CA$eig), 2)

scores_pcoa_2y_bray <- as.data.frame(scores(pcoa_2y_bray, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(betadiv$Bacteria$y2$Meta, by = "Sample") %>%
  mutate(t_m = paste0(Time, "_", Management), Organism = "Bacteria", Dataset = "2-years")

scores_pcoa_2y_bray_16S <- scores_pcoa_2y_bray
variance_2y_bray_16S <- variance_2y_bray[rownames(as.data.frame(variance_2y_bray)) %in% c("MDS1", "MDS2")] %>%
  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("Variance") %>%
  mutate(Variance = "Variance", Organism = "Bacteria", Dataset = "2-years")

p16S_2y <- ggplot(scores_pcoa_2y_bray, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
  theme_pcoa +
  geom_point(shape = 21, stroke = 0.2) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = .1, linewidth = 0.2) +
  labs(x = paste0("PCoA1 ", variance_2y_bray[1], "%"), y = paste0("PCoA2 ", variance_2y_bray[2], "%"),
       title = "Bacteria - 2-years", fill = "Time-Management", color = "Time-Management") +
  scale_color_manual(values = color_timman) +
  scale_fill_manual(values = fill_timman)

ggsave("out/manuscript/betadiv/PCoA_2y_bray_tm_16S.pdf", plot = p16S_2y,
       bg = 'white', width = 2000, height = 1400, dpi = 300, units = "px")

# Bacteria, 3-years -----------------------------------------------------------

pcoa_3y_bray <- capscale(betadiv$Bacteria$y3$Dist$bray ~ 1, distance = "bray")

variance_3y_bray <- round(100 * pcoa_3y_bray$CA$eig / sum(pcoa_3y_bray$CA$eig), 2)

scores_pcoa_3y_bray <- as.data.frame(scores(pcoa_3y_bray, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(betadiv$Bacteria$y3$Meta, by = "Sample") %>%
  mutate(t_m = paste0(Time, "_", Management), Organism = "Bacteria", Dataset = "3-years")

scores_pcoa_3y_bray_16S <- scores_pcoa_3y_bray
variance_3y_bray_16S <- variance_3y_bray[rownames(as.data.frame(variance_3y_bray)) %in% c("MDS1", "MDS2")] %>%
  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("Variance") %>%
  mutate(Variance = "Variance", Organism = "Bacteria", Dataset = "3-years")

p16S_3y <- ggplot(scores_pcoa_3y_bray, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
  theme_pcoa +
  geom_point(shape = 21, stroke = 0.2) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = .1, linewidth = 0.2) +
  labs(x = paste0("PCoA1 ", variance_3y_bray[1], "%"), y = paste0("PCoA2 ", variance_3y_bray[2], "%"),
       title = "Bacteria - 3-years", fill = "Time-Management", color = "Time-Management") +
  scale_color_manual(values = color_timman) +
  scale_fill_manual(values = fill_timman)

ggsave("out/manuscript/betadiv/PCoA_3y_bray_tm_16S.pdf", plot = p16S_3y,
       bg = 'white', width = 2000, height = 1400, dpi = 300, units = "px")

# Fungi, 2-years ----------------------------------------------------------------

pcoa_2y_bray <- capscale(betadiv$Fungi$y2$Dist$bray ~ 1, distance = "bray")

variance_2y_bray <- round(100 * pcoa_2y_bray$CA$eig / sum(pcoa_2y_bray$CA$eig), 2)

scores_pcoa_2y_bray <- as.data.frame(scores(pcoa_2y_bray, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(betadiv$Fungi$y2$Meta, by = "Sample") %>%
  mutate(t_m = paste0(Time, "_", Management), Organism = "Fungi", Dataset = "2-years")

scores_pcoa_2y_bray_ITS <- scores_pcoa_2y_bray
variance_2y_bray_ITS <- variance_2y_bray[rownames(as.data.frame(variance_2y_bray)) %in% c("MDS1", "MDS2")] %>%
  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("Variance") %>%
  mutate(Variance = "Variance", Organism = "Fungi", Dataset = "2-years")

pITS_2y <- ggplot(scores_pcoa_2y_bray, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
  theme_pcoa +
  geom_point(shape = 21, stroke = 0.2) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = .1, linewidth = 0.2) +
  labs(x = paste0("PCoA1 ", variance_2y_bray[1], "%"), y = paste0("PCoA2 ", variance_2y_bray[2], "%"),
       title = "Fungi - 2-years", fill = "Time-Management", color = "Time-Management") +
  scale_color_manual(values = color_timman) +
  scale_fill_manual(values = fill_timman)

ggsave("out/manuscript/betadiv/PCoA_2y_bray_tm_ITS.pdf", plot = pITS_2y,
       bg = 'white', width = 2000, height = 1400, dpi = 300, units = "px")

# Fungi, 3-years ----------------------------------------------------------------

pcoa_3y_bray <- capscale(betadiv$Fungi$y3$Dist$bray ~ 1, distance = "bray")

variance_3y_bray <- round(100 * pcoa_3y_bray$CA$eig / sum(pcoa_3y_bray$CA$eig), 2)

scores_pcoa_3y_bray <- as.data.frame(scores(pcoa_3y_bray, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(betadiv$Fungi$y3$Meta, by = "Sample") %>%
  mutate(t_m = paste0(Time, "_", Management), Organism = "Fungi", Dataset = "3-years")

scores_pcoa_3y_bray_ITS <- scores_pcoa_3y_bray
variance_3y_bray_ITS <- variance_3y_bray[rownames(as.data.frame(variance_3y_bray)) %in% c("MDS1", "MDS2")] %>%
  as.data.frame() %>% t() %>% as.data.frame() %>% rownames_to_column("Variance") %>%
  mutate(Variance = "Variance", Organism = "Fungi", Dataset = "3-years")

pITS_3y <- ggplot(scores_pcoa_3y_bray, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m)) +
  theme_pcoa +
  geom_point(shape = 21, stroke = 0.2) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = .1, linewidth = 0.2) +
  labs(x = paste0("PCoA1 ", variance_3y_bray[1], "%"), y = paste0("PCoA2 ", variance_3y_bray[2], "%"),
       title = "Fungi - 3-years", fill = "Time-Management", color = "Time-Management") +
  scale_color_manual(values = color_timman) +
  scale_fill_manual(values = fill_timman)

ggsave("out/manuscript/betadiv/PCoA_3y_bray_tm_ITS.pdf", plot = pITS_3y,
       bg = 'white', width = 2000, height = 1400, dpi = 300, units = "px")

# Assembled 4-panel manuscript figure (Panel A-D) -------------------------------

(p16S_2y + pITS_2y) / (p16S_3y + pITS_3y) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("out/manuscript/PCoA_bray_panels.pdf",
       bg = 'white', width = 2400, height = 1800, dpi = 240, units = "px")
