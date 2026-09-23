#Script to explore beta-diversity of 2026_AppleMicrobiome ITS data



# Load packages -----------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(readxl)
library(MiscMetabar)
library(vegan)
library(scales)
library(rstatix)



# Define paths ------------------------------------------------------------

path1 <- c("data/")



# Load data ---------------------------------------------------------------

div_obj <- readRDS(paste0(path1, "diversity_data_new.rds"))



# Subset data -------------------------------------------------------------

#whole dataset

meta_all <- div_obj$Fungi$alpha %>%
  mutate(Time = factor(Time, levels = c("May", "July")))


dist_all_bray <- as.matrix(div_obj$Fungi$beta$bray)

dist_all_wuf <- as.matrix(div_obj$Fungi$beta$unifrac_w)


#2-years sampling

meta_2y <- meta_all %>%
  filter(Cultivar %in% c("Honeycrisp", "Spartan"))

dist_2y_bray <- as.dist(dist_all_bray[rownames(dist_all_bray) %in% meta_2y$Sample, colnames(dist_all_bray) %in% meta_2y$Sample])

dist_2y_wuf <- as.dist(dist_all_wuf[rownames(dist_all_wuf) %in% meta_2y$Sample, colnames(dist_all_wuf) %in% meta_2y$Sample])



#2-years, May

meta_2y_may <- meta_2y %>%
  filter(Time == "May")

dist_2y_may_wuf <- as.dist(dist_all_wuf[rownames(dist_all_wuf) %in% meta_2y_may$Sample, colnames(dist_all_wuf) %in% meta_2y_may$Sample])


#2-years, July

meta_2y_july <- meta_2y %>%
  filter(Time == "July")

dist_2y_july_wuf <- as.dist(dist_all_wuf[rownames(dist_all_wuf) %in% meta_2y_july$Sample, colnames(dist_all_wuf) %in% meta_2y_july$Sample])


#3-years sampling

meta_3y <- meta_all %>%
  filter(Cultivar %in% c("Cortland", "Liberty", "Paulared"))

dist_3y_bray <- as.dist(dist_all_bray[rownames(dist_all_bray) %in% meta_3y$Sample, colnames(dist_all_bray) %in% meta_3y$Sample])

dist_3y_wuf <- as.dist(dist_all_wuf[rownames(dist_all_wuf) %in% meta_3y$Sample, colnames(dist_all_wuf) %in% meta_3y$Sample])


#3-years, May 

meta_3y_may <- meta_3y %>%
  filter(Time == "May")

dist_3y_may_wuf <- as.dist(dist_all_wuf[rownames(dist_all_wuf) %in% meta_3y_may$Sample, colnames(dist_all_wuf) %in% meta_3y_may$Sample])


#3-years, July

meta_3y_july <- meta_3y %>%
  filter(Time == "July")

dist_3y_july_wuf <- as.dist(dist_all_wuf[rownames(dist_all_wuf) %in% meta_3y_july$Sample, colnames(dist_all_wuf) %in% meta_3y_july$Sample])


perm_data <- list(meta_2y_may = meta_2y_may, dist_2y_may = dist_2y_may_wuf, meta_2y_july = meta_2y_july, dist_2y_july = dist_2y_july_wuf,
                  meta_3y_may = meta_3y_may, dist_3y_may = dist_3y_may_wuf, meta_3y_july = meta_3y_july, dist_3y_july = dist_3y_july_wuf,
                  meta_2y = meta_2y, dist_2y = dist_2y_wuf, meta_3y = meta_3y, dist_3y = dist_3y_wuf)

saveRDS(perm_data, "Permanova_input_ITS.rds")



# PCoA and PERMANOVA whole dataset (Bray-Curtis) ----------------------------------------

set.seed(230726)


pcoa_all_bray <- capscale(dist_all_bray~1, distance = "bray")

variance_all_bray <- round(100 * pcoa_all_bray$CA$eig / sum(pcoa_all_bray$CA$eig), 2)

scores_pcoa_all_bray <- as.data.frame(scores(pcoa_all_bray, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(meta_all, by = "Sample") %>%
  mutate(y_t = paste0(Year, "_", Time), Year = as.character(Year), t_m = paste0(Time, "_", Management)) 


ggplot(scores_pcoa_all_bray, aes(x = MDS1, y = MDS2, color = Management, fill = Management, shape = Management)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Management)) +
  labs(x = paste0("PCoA1 ", variance_all_bray[1], "%"), y = paste0("PCoA2 ", variance_all_bray[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_man, values = colors_man) +
  scale_fill_manual(breaks = labels_man, values = fill_man) +
  scale_shape_manual(breaks = labels_man, values = shapes_man)

ggsave("Fig_bray/PCoA_all_bray_time_man_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_bray, aes(x = MDS1, y = MDS2, color = Management, fill = Management, shape = Time)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Management)) +
  labs(x = paste0("PCoA1 ", variance_all_bray[1], "%"), y = paste0("PCoA2 ", variance_all_bray[2], "%")) +
  scale_color_manual(breaks = labels_man, values = colors_man) +
  scale_fill_manual(breaks = labels_man, values = fill_man) +
  scale_shape_manual(breaks = labels_time, values = shapes_time2)

ggsave("Fig_bray/PCoA_all_bray_man_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_bray, aes(x = MDS1, y = MDS2, color = Cultivar, fill = Cultivar, shape = Cultivar)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Cultivar)) +
  labs(x = paste0("PCoA1 ", variance_all_bray[1], "%"), y = paste0("PCoA2 ", variance_all_bray[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_cult, values = colors_cult) +
  scale_fill_manual(breaks = labels_cult, values = fill_cult) +
  scale_shape_manual(breaks = labels_cult, values = shapes_cult)

ggsave("Fig_bray/PCoA_all_bray_time_cult_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_bray, aes(x = MDS1, y = MDS2, color = Location, fill = Location, shape = Location)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Location)) +
  labs(x = paste0("PCoA1 ", variance_all_bray[1], "%"), y = paste0("PCoA2 ", variance_all_bray[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_loc, values = colors_loc) +
  scale_fill_manual(breaks = labels_loc, values = fill_loc) +
  scale_shape_manual(breaks = labels_loc, values = shapes_loc)

ggsave("Fig_bray/PCoA_all_bray_time_loc_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_bray, aes(x = MDS1, y = MDS2, color = Year, fill = Year, shape = Year)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Year)) +
  labs(x = paste0("PCoA1 ", variance_all_bray[1], "%"), y = paste0("PCoA2 ", variance_all_bray[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_year, values = colors_year) +
  scale_fill_manual(breaks = labels_year, values = fill_year) +
  scale_shape_manual(breaks = labels_year, values = shapes_year)

ggsave("Fig_bray/PCoA_all_bray_time_year_ITS.png", dpi = 300)



# PCoA and PERMANOVA whole dataset (weighted UniFrac) ---------------------

set.seed(230726)


pcoa_all_wuf <- capscale(dist_all_wuf~1, distance = "wunifrac")

variance_all_wuf <- round(100 * pcoa_all_wuf$CA$eig / sum(pcoa_all_wuf$CA$eig), 2)

scores_pcoa_all_wuf <- as.data.frame(scores(pcoa_all_wuf, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(meta_all, by = "Sample") %>%
  mutate(y_t = paste0(Year, "_", Time), Year = as.character(Year), t_m = paste0(Time, "_", Management)) 


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, color = Management, fill = Management, shape = Management)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Management)) +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_man, values = colors_man) +
  scale_fill_manual(breaks = labels_man, values = fill_man) +
  scale_shape_manual(breaks = labels_man, values = shapes_man)

ggsave("Fig_wuf/PCoA_all_wuf_time_man_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m, shape = t_m)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = t_m)) +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  scale_color_manual(breaks = labels_timman, values = colors_timman) +
  scale_fill_manual(breaks = labels_timman, values = fill_timman) +
  scale_shape_manual(breaks = labels_timman, values = shapes_timman)

ggsave("Fig_wuf/PCoA_all_wuf_tm_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, color = Management, fill = Management, shape = Time)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Management)) +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  scale_color_manual(breaks = labels_man, values = colors_man) +
  scale_fill_manual(breaks = labels_man, values = fill_man) +
  scale_shape_manual(breaks = labels_time, values = shapes_time2)

ggsave("Fig_wuf/PCoA_all_wuf_man_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, color = Cultivar, fill = Cultivar, shape = Cultivar)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Cultivar)) +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_cult, values = colors_cult) +
  scale_fill_manual(breaks = labels_cult, values = fill_cult) +
  scale_shape_manual(breaks = labels_cult, values = shapes_cult)

ggsave("Fig_wuf/PCoA_all_wuf_time_cult_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, color = Location, fill = Location, shape = Location)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Location)) +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_loc, values = colors_loc) +
  scale_fill_manual(breaks = labels_loc, values = fill_loc) +
  scale_shape_manual(breaks = labels_loc, values = shapes_loc)

ggsave("Fig_wuf/PCoA_all_wuf_time_loc_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, color = Year, fill = Year, shape = Year)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Year)) +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_year, values = colors_year) +
  scale_fill_manual(breaks = labels_year, values = fill_year) +
  scale_shape_manual(breaks = labels_year, values = shapes_year)

ggsave("Fig_wuf/PCoA_all_wuf_time_year_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, color = Time, fill = Time, shape = Time)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Time)) +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  scale_color_manual(breaks = labels_time, values = colors_time) +
  scale_fill_manual(breaks = labels_time, values = fill_time) +
  scale_shape_manual(breaks = labels_time, values = shapes_time)

ggsave("Fig_wuf/PCoA_all_wuf_time_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, fill = mean_temp)) +
  theme_pcoa +
  geom_point(size = 5, shape = 21, color = "black") +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  scale_fill_gradient(high = grad_low, low = grad_hi, limits = c(11, 24))

ggsave("Fig_wuf/PCoA_all_wuf_temp_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, fill = precip_72h)) +
  theme_pcoa +
  geom_point(size = 5, shape = 21, color = "black") +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  scale_fill_gradient(high = grad_low, low = grad_hi, limits = c(1, 86))

ggsave("Fig_wuf/PCoA_all_wuf_precip_ITS.png", dpi = 300)


ggplot(scores_pcoa_all_wuf, aes(x = MDS1, y = MDS2, fill = deg_h)) +
  theme_pcoa +
  geom_point(size = 5, shape = 21, color = "black") +
  labs(x = paste0("PCoA1 ", variance_all_wuf[1], "%"), y = paste0("PCoA2 ", variance_all_wuf[2], "%")) +
  scale_fill_gradient(high = grad_low, low = grad_hi, limits = c(228, 2286))

ggsave("Fig_wuf/PCoA_all_wuf_degh_ITS.png", dpi = 300)



# PCoA and PERMANOVA 2-Years (Bray-Curtis) --------------------------------

set.seed(230726)


pcoa_2y_bray <- capscale(dist_2y_bray~1, distance = "bray")

variance_2y_bray <- round(100 * pcoa_2y_bray$CA$eig / sum(pcoa_2y_bray$CA$eig), 2)

scores_pcoa_2y_bray <- as.data.frame(scores(pcoa_2y_bray, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(meta_2y, by = "Sample") %>%
  mutate(y_t = paste0(Year, "_", Time), Year = as.character(Year), t_m = paste0(Time, "_", Management)) 


ggplot(scores_pcoa_2y_bray, aes(x = MDS1, y = MDS2, color = Management, fill = Management, shape = Management)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Management)) +
  labs(x = paste0("PCoA1 ", variance_2y_bray[1], "%"), y = paste0("PCoA2 ", variance_2y_bray[2], "%")) +
  facet_wrap(~Time)  +
  scale_color_manual(breaks = labels_man, values = colors_man) +
  scale_fill_manual(breaks = labels_man, values = fill_man) +
  scale_shape_manual(breaks = labels_man, values = shapes_man)

ggsave("Fig_bray/PCoA_2y_time_man_ITS.png", dpi = 300)


ggplot(scores_pcoa_2y_bray, aes(x = MDS1, y = MDS2, color = Cultivar, fill = Cultivar, shape = Cultivar)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Cultivar)) +
  labs(x = paste0("PCoA1 ", variance_2y_bray[1], "%"), y = paste0("PCoA2 ", variance_2y_bray[2], "%")) +
  facet_wrap(~Time)  +
  scale_color_manual(breaks = labels_cult, values = colors_cult) +
  scale_fill_manual(breaks = labels_cult, values = fill_cult) +
  scale_shape_manual(breaks = labels_cult, values = shapes_cult)

ggsave("Fig_bray/PCoA_2y_time_cult_ITS.png", dpi = 300)


ggplot(scores_pcoa_2y_bray, aes(x = MDS1, y = MDS2, color = Location, fill = Location, shape = Location)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Location)) +
  labs(x = paste0("PCoA1 ", variance_2y_bray[1], "%"), y = paste0("PCoA2 ", variance_2y_bray[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_loc, values = colors_loc) +
  scale_fill_manual(breaks = labels_loc, values = fill_loc) +
  scale_shape_manual(breaks = labels_loc, values = shapes_loc)

ggsave("Fig_bray/PCoA_2y_time_loc_ITS.png", dpi = 300)



disp_loc_2y_bray <- betadisper(dist_2y_bray, meta_2y$Location, type ="centroid")
anova_loc_2y_bray <- as.data.frame(anova(disp_loc_2y_bray))

disp_man_2y_bray <- betadisper(dist_2y_bray, meta_2y$Management, type ="centroid")
anova_man_2y_bray <- as.data.frame(anova(disp_man_2y_bray))

disp_year_2y_bray <- betadisper(dist_2y_bray, meta_2y$Year, type ="centroid")
anova_year_2y_bray <- as.data.frame(anova(disp_year_2y_bray))

disp_cult_2y_bray <- betadisper(dist_2y_bray, meta_2y$Cultivar, type ="centroid")
anova_cult_2y_bray <- as.data.frame(anova(disp_cult_2y_bray))



perm_loc_man_cult_2y_bray <- adonis2(dist_2y_bray~Management*Location*Time+Cultivar+Seq_depth, meta_2y, permutations = 999, method = "bray", by = "terms")

perm_loc_man_cult_2y_bray %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_bray/Perm_2y_loc_man_time_ITS.csv")

perm_temp_prec_dh_man_cul_2y_bray <- adonis2(dist_2y_bray~mean_temp+precip_72h+deg_h+Management+Cultivar+Seq_depth, meta_2y, permutations = 999, method = "bray", by = "terms")

perm_temp_prec_dh_man_cul_2y_bray %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_bray/Perm_2y_temp_prec_dh_man_ITS.csv")



# PCoA and PERMANOVA 2-Years (weighted UniFrac) ---------------------------

set.seed(230726)


pcoa_2y_wuf <- capscale(dist_2y_wuf~1, distance = "wunifrac")

variance_2y_wuf <- round(100 * pcoa_2y_wuf$CA$eig / sum(pcoa_2y_wuf$CA$eig), 2)

variance_2y_wuf_ITS <- as.data.frame(variance_2y_wuf[rownames(as.data.frame(variance_2y_wuf)) %in% c("MDS1", "MDS2")]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Variance") %>%
  mutate(Variance = paste0("Variance"), Organism = "Fungi", Dataset = "2-years")


scores_pcoa_2y_wuf <- as.data.frame(scores(pcoa_2y_wuf, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(meta_2y, by = "Sample") %>%
  mutate(y_t = paste0(Year, "_", Time), Year = as.character(Year), t_m = paste0(Time, "_", Management),
         Organism = paste0("Fungi"), Dataset = paste0("2-years")) 

scores_pcoa_2y_wuf_ITS <- scores_pcoa_2y_wuf


ggplot(scores_pcoa_2y_wuf, aes(x = MDS1, y = MDS2, color = Management, fill = Management, shape = Management)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Management)) +
  labs(x = paste0("PCoA1 ", variance_2y_wuf[1], "%"), y = paste0("PCoA2 ", variance_2y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_man, values = colors_man) +
  scale_fill_manual(breaks = labels_man, values = fill_man) +
  scale_shape_manual(breaks = labels_man, values = shapes_man)

ggsave("Fig_wuf/PCoA_2y_time_man_ITS.png", dpi = 300)


pITS_2y <- ggplot(scores_pcoa_2y_wuf, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m, shape = t_m)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = t_m)) +
  labs(x = paste0("PCoA1 ", variance_2y_wuf[1], "%"), y = paste0("PCoA2 ", variance_2y_wuf[2], "%"), title = "Fungi - 2-years", fill = "Time-Management", shape = "Time-Management", color = "Time-Management")  +
  scale_color_manual(breaks = labels_timman, values = colors_timman) +
  scale_fill_manual(breaks = labels_timman, values = fill_timman) +
  scale_shape_manual(breaks = labels_timman, values = shapes_timman)

ggsave("Fig_wuf/PCoA_2y_wuf_tm_ITS.png", dpi = 300)


ggplot(scores_pcoa_2y_wuf, aes(x = MDS1, y = MDS2, color = Cultivar, fill = Cultivar, shape = Cultivar)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Cultivar)) +
  labs(x = paste0("PCoA1 ", variance_2y_wuf[1], "%"), y = paste0("PCoA2 ", variance_2y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_cult, values = colors_cult) +
  scale_fill_manual(breaks = labels_cult, values = fill_cult) +
  scale_shape_manual(breaks = labels_cult, values = shapes_cult)

ggsave("Fig_wuf/PCoA_2y_time_cult_ITS.png", dpi = 300)


ggplot(scores_pcoa_2y_wuf, aes(x = MDS1, y = MDS2, color = Location, fill = Location, shape = Location)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Location)) +
  labs(x = paste0("PCoA1 ", variance_2y_wuf[1], "%"), y = paste0("PCoA2 ", variance_2y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_loc, values = colors_loc) +
  scale_fill_manual(breaks = labels_loc, values = fill_loc) +
  scale_shape_manual(breaks = labels_loc, values = shapes_loc)

ggsave("Fig_wuf/PCoA_2y_time_loc_ITS.png", dpi = 300)


ggplot(scores_pcoa_2y_wuf, aes(x = MDS1, y = MDS2, fill = mean_temp)) +
  theme_pcoa +
  geom_point(size = 5, shape = 21, color = "black") +
  labs(x = paste0("PCoA1 ", variance_2y_wuf[1], "%"), y = paste0("PCoA2 ", variance_2y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_fill_gradient(high = grad_low, low = grad_hi, limits = c(12, 24))

ggsave("Fig_wuf/PCoA_2y_time_temp_ITS.png", dpi = 300)


ggplot(scores_pcoa_2y_wuf, aes(x = MDS1, y = MDS2, fill = precip_72h)) +
  theme_pcoa +
  geom_point(size = 5, shape = 21, color = "black") +
  labs(x = paste0("PCoA1 ", variance_2y_wuf[1], "%"), y = paste0("PCoA2 ", variance_2y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_fill_gradient(high = grad_low, low = grad_hi, limits = c(1, 86))

ggsave("Fig_wuf/PCoA_2y_time_prec_ITS.png", dpi = 300)



disp_loc_2y_wuf <- betadisper(dist_2y_wuf, meta_2y$Location, type ="centroid")
anova_loc_2y_wuf <- as.data.frame(anova(disp_loc_2y_wuf)) #0.000138098

disp_man_2y_wuf <- betadisper(dist_2y_wuf, meta_2y$Management, type ="centroid")
anova_man_2y_wuf <- as.data.frame(anova(disp_man_2y_wuf)) #ns

disp_year_2y_wuf <- betadisper(dist_2y_wuf, meta_2y$Year, type ="centroid")
anova_year_2y_wuf <- as.data.frame(anova(disp_year_2y_wuf)) #0.01307386

disp_cult_2y_wuf <- betadisper(dist_2y_wuf, meta_2y$Cultivar, type ="centroid")
anova_cult_2y_wuf <- as.data.frame(anova(disp_cult_2y_wuf)) #ns

disp_time_2y_wuf <- betadisper(dist_2y_wuf, meta_2y$Time, type ="centroid")
anova_time_2y_wuf <- as.data.frame(anova(disp_time_2y_wuf)) #0.0143386

disp_temp_2y_wuf <- betadisper(dist_2y_wuf, meta_2y$mean_temp, type ="centroid")
anova_temp_2y_wuf <- as.data.frame(anova(disp_temp_2y_wuf)) #6.198143e-06

disp_precip_2y_wuf <- betadisper(dist_2y_wuf, meta_2y$precip_72h_z, type ="centroid")
anova_precip_2y_wuf <- as.data.frame(anova(disp_precip_2y_wuf)) #6.198143e-06

disp_degh_2y_wuf <- betadisper(dist_2y_wuf, meta_2y$deg_h_z, type ="centroid")
anova_degh_2y_wuf <- as.data.frame(anova(disp_degh_2y_wuf)) #6.198143e-06



perm_loc_man_cult_2y_wuf <- adonis2(dist_2y_wuf~seq_depth_z+Management*Location*Time+Cultivar, meta_2y, permutations = 999, method = "wunifrac", by = "terms")

perm_loc_man_cult_2y_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_2y_loc_man_time_ITS.csv")

perm_loc_man_cult_2y_may_wuf <- adonis2(dist_2y_may_wuf~seq_depth_z+Management*Location+Cultivar, meta_2y_may, permutations = 999, method = "wunifrac", by = "terms")

perm_loc_man_cult_2y_may_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_2y_may_loc_man_time_ITS.csv")


perm_loc_man_cult_2y_july_wuf <- adonis2(dist_2y_july_wuf~seq_depth_z+Management*Location+Cultivar, meta_2y_july, permutations = 999, method = "wunifrac", by = "terms")

perm_loc_man_cult_2y_july_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_2y_july_loc_man_time_ITS.csv")


perm_temp_prec_dh_man_cul_2y_wuf <- adonis2(dist_2y_wuf~seq_depth_z+mean_temp+precip_72h_z+deg_h_z+Management+Cultivar, meta_2y, permutations = 999, method = "wunifrac", by = "terms")

perm_temp_prec_dh_man_cul_2y_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_2y_temp_prec_dh_man_ITS.csv")

perm_temp_prec_dh_man_cul_2y_may_wuf <- adonis2(dist_2y_may_wuf~seq_depth_z+mean_temp+precip_72h_z+deg_h_z+Management+Cultivar, meta_2y_may, permutations = 999, method = "wunifrac", by = "terms")

perm_temp_prec_dh_man_cul_2y_may_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_2y_may_temp_prec_dh_man_ITS.csv")


perm_temp_prec_dh_man_cul_2y_july_wuf <- adonis2(dist_2y_july_wuf~seq_depth_z+mean_temp+precip_72h_z+deg_h_z+Management+Cultivar, meta_2y_july, permutations = 999, method = "wunifrac", by = "terms")

perm_temp_prec_dh_man_cul_2y_july_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_2y_july_temp_prec_dh_man_ITS.csv")


# PCoA and PERMANOVA 3-Years (Bray-Curtis) --------------------------------

set.seed(230726)


pcoa_3y_bray <- capscale(dist_3y_bray~1, distance = "bray")

variance_3y_bray <- round(100 * pcoa_3y_bray$CA$eig / sum(pcoa_3y_bray$CA$eig), 2)

scores_pcoa_3y_bray <- as.data.frame(scores(pcoa_3y_bray, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(meta_3y, by = "Sample") %>%
  mutate(y_t = paste0(Year, "_", Time), Year = as.character(Year), t_m = paste0(Time, "_", Management)) 


ggplot(scores_pcoa_3y_bray, aes(x = MDS1, y = MDS2, color = Management, fill = Management, shape = Management)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Management)) +
  labs(x = paste0("PCoA1 ", variance_3y_bray[1], "%"), y = paste0("PCoA2 ", variance_3y_bray[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_man, values = colors_man) +
  scale_fill_manual(breaks = labels_man, values = fill_man) +
  scale_shape_manual(breaks = labels_man, values = shapes_man)

ggsave("Fig_bray/PCoA_3y_time_man_ITS.png", dpi = 300)


ggplot(scores_pcoa_3y_bray, aes(x = MDS1, y = MDS2, color = Cultivar, fill = Cultivar, shape = Cultivar)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Cultivar)) +
  labs(x = paste0("PCoA1 ", variance_3y_bray[1], "%"), y = paste0("PCoA2 ", variance_3y_bray[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_cult, values = colors_cult) +
  scale_fill_manual(breaks = labels_cult, values = fill_cult) +
  scale_shape_manual(breaks = labels_cult, values = shapes_cult)

ggsave("Fig_bray/PCoA_3y_time_cult_ITS.png", dpi = 300)


ggplot(scores_pcoa_3y_bray, aes(x = MDS1, y = MDS2, color = Location, fill = Location, shape = Location)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Location)) +
  labs(x = paste0("PCoA1 ", variance_3y_bray[1], "%"), y = paste0("PCoA2 ", variance_3y_bray[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_loc, values = colors_loc) +
  scale_fill_manual(breaks = labels_loc, values = fill_loc) +
  scale_shape_manual(breaks = labels_loc, values = shapes_loc)

ggsave("Fig_bray/PCoA_3y_time_loc_ITS.png", dpi = 300)



disp_loc_3y_bray <- betadisper(dist_3y_bray, meta_3y$Location, type ="centroid")
anova_loc_3y_bray <- as.data.frame(anova(disp_loc_3y_bray))

disp_man_3y_bray <- betadisper(dist_3y_bray, meta_3y$Management, type ="centroid")
anova_man_3y_bray <- as.data.frame(anova(disp_man_3y_bray))

disp_Year_3y_bray <- betadisper(dist_3y_bray, meta_3y$Year, type ="centroid")
anova_Year_3y_bray <- as.data.frame(anova(disp_Year_3y_bray))

disp_cult_3y_bray <- betadisper(dist_3y_bray, meta_3y$Cultivar, type ="centroid")
anova_cult_3y_bray <- as.data.frame(anova(disp_cult_3y_bray))


perm_loc_man_cult_3y_bray <- adonis2(dist_3y_bray~Management*Location*Time+Cultivar+Seq_depth, meta_3y, permutations = 999, method = "wunifrac", by = "terms")

perm_loc_man_cult_3y_bray %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_bray/Perm_3y_loc_man_time_ITS.csv")

perm_temp_prec_dh_man_cul_3y_bray <- adonis2(dist_3y_bray~mean_temp+precip_72h+deg_h+Management+Cultivar+Seq_depth, meta_3y, permutations = 999, method = "wunifrac", by = "terms")

perm_temp_prec_dh_man_cul_3y_bray %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_bray/Perm_3y_temp_prec_dh_man_ITS.csv")



# PCoA and PERMANOVA 3-Years (weighted UniFrac) ---------------------------

set.seed(230726)


pcoa_3y_wuf <- capscale(dist_3y_wuf~1, distance = "wunifrac")

variance_3y_wuf <- round(100 * pcoa_3y_wuf$CA$eig / sum(pcoa_3y_wuf$CA$eig), 2)

variance_3y_wuf_ITS <- as.data.frame(variance_3y_wuf[rownames(as.data.frame(variance_3y_wuf)) %in% c("MDS1", "MDS2")]) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Variance") %>%
  mutate(Variance = paste0("Variance"), Organism = "Fungi", Dataset = "3-years")


scores_pcoa_3y_wuf <- as.data.frame(scores(pcoa_3y_wuf, display = "sites")) %>%
  rownames_to_column("Sample") %>%
  merge(meta_3y, by = "Sample") %>%
  mutate(y_t = paste0(Year, "_", Time), Year = as.character(Year), t_m = paste0(Time, "_", Management),
         Organism = paste0("Fungi"), Dataset = paste0("3-years"))

scores_pcoa_3y_wuf_ITS <- scores_pcoa_3y_wuf


ggplot(scores_pcoa_3y_wuf, aes(x = MDS1, y = MDS2, color = Management, fill = Management, shape = Management)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Management)) +
  labs(x = paste0("PCoA1 ", variance_3y_wuf[1], "%"), y = paste0("PCoA2 ", variance_3y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_man, values = colors_man) +
  scale_fill_manual(breaks = labels_man, values = fill_man) +
  scale_shape_manual(breaks = labels_man, values = shapes_man)

ggsave("Fig_wuf/PCoA_3y_time_man_ITS.png", dpi = 300)


pITS_3y <- ggplot(scores_pcoa_3y_wuf, aes(x = MDS1, y = MDS2, color = t_m, fill = t_m, shape = t_m)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = t_m)) +
  labs(x = paste0("PCoA1 ", variance_3y_wuf[1], "%"), y = paste0("PCoA2 ", variance_3y_wuf[2], "%"), title = "Fungi - 3-years", fill = "Time-Management", shape = "Time-Management", color = "Time-Management") +
  scale_color_manual(breaks = labels_timman, values = colors_timman) +
  scale_fill_manual(breaks = labels_timman, values = fill_timman) +
  scale_shape_manual(breaks = labels_timman, values = shapes_timman)

ggsave("Fig_wuf/PCoA_3y_wuf_tm_ITS.png", dpi = 300)


ggplot(scores_pcoa_3y_wuf, aes(x = MDS1, y = MDS2, color = Cultivar, fill = Cultivar, shape = Cultivar)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Cultivar)) +
  labs(x = paste0("PCoA1 ", variance_3y_wuf[1], "%"), y = paste0("PCoA2 ", variance_3y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_cult, values = colors_cult) +
  scale_fill_manual(breaks = labels_cult, values = fill_cult) +
  scale_shape_manual(breaks = labels_cult, values = shapes_cult)

ggsave("Fig_wuf/PCoA_3y_time_cult_ITS.png", dpi = 300)


ggplot(scores_pcoa_3y_wuf, aes(x = MDS1, y = MDS2, color = Location, fill = Location, shape = Location)) +
  theme_pcoa +
  geom_point(size = 5) +
  stat_ellipse(level=0.95, geom = "polygon", alpha = .1, aes(group = Location)) +
  labs(x = paste0("PCoA1 ", variance_3y_wuf[1], "%"), y = paste0("PCoA2 ", variance_3y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_color_manual(breaks = labels_loc, values = colors_loc) +
  scale_fill_manual(breaks = labels_loc, values = fill_loc) +
  scale_shape_manual(breaks = labels_loc, values = shapes_loc)

ggsave("Fig_wuf/PCoA_3y_time_loc_ITS.png", dpi = 300)


ggplot(scores_pcoa_3y_wuf, aes(x = MDS1, y = MDS2, fill = mean_temp)) +
  theme_pcoa +
  geom_point(size = 5, shape = 21, color = "black") +
  labs(x = paste0("PCoA1 ", variance_3y_wuf[1], "%"), y = paste0("PCoA2 ", variance_3y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_fill_gradient(high = grad_low, low = grad_hi, limits = c(11, 24))

ggsave("Fig_wuf/PCoA_3y_time_temp_ITS.png", dpi = 300)


ggplot(scores_pcoa_3y_wuf, aes(x = MDS1, y = MDS2, fill = precip_72h)) +
  theme_pcoa +
  geom_point(size = 5, shape = 21, color = "black") +
  labs(x = paste0("PCoA1 ", variance_3y_wuf[1], "%"), y = paste0("PCoA2 ", variance_3y_wuf[2], "%")) +
  facet_wrap(~Time) +
  scale_fill_gradient(high = grad_low, low = grad_hi, limits = c(1, 86))

ggsave("Fig_wuf/PCoA_3y_time_prec_ITS.png", dpi = 300)



disp_loc_3y_wuf <- betadisper(dist_3y_wuf, meta_3y$Location, type ="centroid")
anova_loc_3y_wuf <- as.data.frame(anova(disp_loc_3y_wuf)) #6.031411e-06

disp_man_3y_wuf <- betadisper(dist_3y_wuf, meta_3y$Management, type ="centroid")
anova_man_3y_wuf <- as.data.frame(anova(disp_man_3y_wuf)) #2.641732e-07

disp_Year_3y_wuf <- betadisper(dist_3y_wuf, meta_3y$Year, type ="centroid")
anova_Year_3y_wuf <- as.data.frame(anova(disp_Year_3y_wuf)) #ns

disp_cult_3y_wuf <- betadisper(dist_3y_wuf, meta_3y$Cultivar, type ="centroid")
anova_cult_3y_wuf <- as.data.frame(anova(disp_cult_3y_wuf)) #ns

disp_time_3y_wuf <- betadisper(dist_3y_wuf, meta_3y$Time, type ="centroid")
anova_time_3y_wuf <- as.data.frame(anova(disp_time_3y_wuf)) #0.009001633

disp_temp_3y_wuf <- betadisper(dist_3y_wuf, meta_3y$mean_temp, type ="centroid")
anova_temp_3y_wuf <- as.data.frame(anova(disp_temp_3y_wuf)) #4.321043e-10

disp_precip_3y_wuf <- betadisper(dist_3y_wuf, meta_3y$precip_72h_z, type ="centroid")
anova_precip_3y_wuf <- as.data.frame(anova(disp_precip_3y_wuf)) #3.652386e-11

disp_degh_3y_wuf <- betadisper(dist_3y_wuf, meta_3y$deg_h_z, type ="centroid")
anova_degh_3y_wuf <- as.data.frame(anova(disp_degh_3y_wuf)) #4.321043e-10


perm_loc_man_cult_3y_wuf <- adonis2(dist_3y_wuf~seq_depth_z+Management*Location*Time+Cultivar, meta_3y, permutations = 999, method = "wunifrac", by = "terms")

perm_loc_man_cult_3y_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_3y_loc_man_time_ITS.csv")

perm_loc_man_cult_3y_may_wuf <- adonis2(dist_3y_may_wuf~seq_depth_z+Management*Location+Cultivar, meta_3y_may, permutations = 999, method = "wunifrac", by = "terms")

perm_loc_man_cult_3y_may_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_3y_may_loc_man_time_ITS.csv")


perm_loc_man_cult_3y_july_wuf <- adonis2(dist_3y_july_wuf~seq_depth_z+Management*Location+Cultivar, meta_3y_july, permutations = 999, method = "wunifrac", by = "terms")

perm_loc_man_cult_3y_july_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_3y_july_loc_man_time_ITS.csv")


perm_temp_prec_dh_man_cul_3y_wuf <- adonis2(dist_3y_wuf~seq_depth_z+mean_temp+precip_72h_z+deg_h_z+Management+Cultivar, meta_3y, permutations = 999, method = "wunifrac", by = "terms")

perm_temp_prec_dh_man_cul_3y_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_3y_temp_prec_dh_man_ITS.csv")

perm_temp_prec_dh_man_cul_3y_may_wuf <- adonis2(dist_3y_may_wuf~seq_depth_z+mean_temp+precip_72h_z+deg_h_z+Management+Cultivar, meta_3y_may, permutations = 999, method = "wunifrac", by = "terms")

perm_temp_prec_dh_man_cul_3y_may_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_3y_may_temp_prec_dh_man_ITS.csv")


perm_temp_prec_dh_man_cul_3y_july_wuf <- adonis2(dist_3y_july_wuf~seq_depth_z+mean_temp+precip_72h_z+deg_h_z+Management+Cultivar, meta_3y_july, permutations = 999, method = "wunifrac", by = "terms")

perm_temp_prec_dh_man_cul_3y_july_wuf %>%
  rownames_to_column("Parameter") %>%
  write_csv("Perm_wuf/Perm_3y_july_temp_prec_dh_man_ITS.csv")




