# Beta-diversity statistics (weighted UniFrac), both barcodes
# - PERMANOVA (adonis2) on all data and per subset x month
# - dispersion homogeneity (betadisper) for each PERMANOVA factor

pacman::p_load(tidyverse, vegan, update = FALSE)
set.seed(230726)

#source('src/0.0_Config.R')
source('src/utils/beta_div_saver.R') # save_stat_xlsx()
source('src/0.0_ggplot_themes.R')

stats_dir_out <- 'out/manuscript'
betadiv <- readRDS('data/diversity_subsets.rds')

perms <- list()

# BACTERIA ---- 

perms[["All data 16S"]] <- adonis2(
  betadiv$Bacteria$all$Dist$wuf ~ seq_depth_z + Year*Time*Location,
  betadiv$Bacteria$all$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["3Y - May 16S"]] <- adonis2(
  betadiv$Bacteria$y3_may$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Location, 
  betadiv$Bacteria$y3_may$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["3Y - July 16S"]] <- adonis2(
  betadiv$Bacteria$y3_july$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Location, 
  betadiv$Bacteria$y3_july$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["2Y - May 16S"]] <- adonis2(
  betadiv$Bacteria$y2_may$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Location * Management, 
  betadiv$Bacteria$y2_may$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["2Y - July 16S"]] <- adonis2(
  betadiv$Bacteria$y2_july$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Location * Management, 
  betadiv$Bacteria$y2_july$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")


# FUNGI -----------

perms[["All data ITS"]]<- adonis2(
  betadiv$Fungi$all$Dist$wuf ~ seq_depth_z + Year*Time*Location,
  betadiv$Fungi$all$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["3Y - May ITS"]] <- adonis2(
  betadiv$Fungi$y3_may$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Location, 
  betadiv$Fungi$y3_may$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["3Y - July ITS"]] <- adonis2(
  betadiv$Fungi$y3_july$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z  + Year * Location, 
  betadiv$Fungi$y3_july$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["2Y - May ITS"]] <- adonis2(
  betadiv$Fungi$y2_may$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Location * Management, 
  betadiv$Fungi$y2_may$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["2Y - July ITS"]] <- adonis2(
  betadiv$Fungi$y2_july$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Location * Management, 
  betadiv$Fungi$y2_july$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")


# Export: one sheet per perms element ------
save_stat_xlsx("out/manuscript/betadiv_permanova.xlsx", perms)

# DISPERSION (betadisper) ------------------
# Homogeneity of dispersion for each categorical PERMANOVA factor, same subsets
# - a significant PERMANOVA term can reflect spread rather than centroid shift
# - centroid-based (as PERMANOVA), permutation test
# - continuous covariates not tested (betadisper needs groups)

disp_models <- tribble(
  ~Model,      ~slot,     ~factors,
  "All data",  "all",     c("Year", "Time", "Location"),
  "3Y - May",  "y3_may",  c("Year", "Location"),
  "3Y - July", "y3_july", c("Year", "Location"),
  "2Y - May",  "y2_may",  c("Year", "Location", "Management"),
  "2Y - July", "y2_july", c("Year", "Location", "Management")
)

disp_specs <- disp_models %>%
  crossing(Kingdom = c("Bacteria", "Fungi")) %>%
  unnest_longer(factors, values_to = "Factor")

set.seed(230726); dispersion <- disp_specs %>%
  mutate(res = pmap(list(Kingdom, slot, Factor), \(k, s, f) {
    sub <- betadiv[[k]][[s]]
    stopifnot(identical(labels(sub$Dist$wuf), as.character(sub$Meta$Sample)))
    grp <- droplevels(factor(sub$Meta[[f]]))
    bd  <- betadisper(sub$Dist$wuf, grp, type = "centroid")
    pt  <- permutest(bd, permutations = 999, parallel = 6)
    tibble(
      F = pt$tab$F[1],
      p = pt$tab$`Pr(>F)`[1],
      # mean distance to group centroid: which group is more dispersed
      `Mean distance to centroid` = paste(
        sprintf("%s: %.2f", levels(grp), tapply(bd$distances, grp, mean)), collapse = "; "))
  })) %>%
  unnest(res) %>%
  mutate(`p (BH)` = round(p.adjust(p, method = "BH"), 3),
         Kingdom = factor(Kingdom, levels = c("Bacteria", "Fungi")),
         Model   = factor(Model, levels = disp_models$Model)) %>%
  arrange(Kingdom, Model) %>%
  select(Kingdom, Model, Factor, `p (BH)`, `Mean distance to centroid`)

save_stat_kable("out/manuscript/betadiv_betadisper.html", dispersion, rownames_to = NULL,
  caption = paste(
    "Homogeneity of multivariate dispersion (betadisper, weighted UniFrac, distance to",
    "group centroid; permutation test, 999 permutations) for each factor of the",
    "corresponding PERMANOVA model. p (BH): Benjamini-Hochberg across all tests."))

# A PLOT ------------------

model_vars <- c(
  'Residual',
  'seq_depth_z', 'deg_h_z', 'precip_72h_z', 'mean_temp', 
  'Year', 'Location', 'Management', 'Time',
  'Year:Location', 'Year:Location:Management', 'Year:Management', 
  'Location:Management',
  'Year:Time', 'Year:Time:Location', 'Time:Location'
)

perm_out_full <- imap(perms, function(permanova, model_name){
  res <- perms[[model_name]]
  tibble(
    Model = model_name,
    variable = rownames(res),
    R2 = res$R2,
    p = res$`Pr(>F)`,
    df = res$Df,
    pseudoF = res$F
  ) 
}) %>% list_rbind() %>% 
  filter(variable != 'Total') %>% 
  mutate(
    variable = factor(
      variable, levels = model_vars
    ),
    Kingdom = case_when(
      str_detect(Model, "16S") ~ 'Bacteria', TRUE ~ 'Fungi'
    ),
    Model = str_remove(Model, "\\s+\\S*$"),
    
  )
# 
# perm_out_full %<>% 
#   mutate(variable = factor(
#     variable, levels = c(model_vars, 'Residual')))

palette <- c('grey90', MetBrewer::met.brewer('Redon', n=length(model_vars)-1))
perm_out_full %>% 
  ggplot(aes(y = Model, x = R2, #alpha = p_sig,
             fill = variable)) +
  geom_col() +
  scale_fill_manual(values = palette) +
  facet_grid(Kingdom~., scale = 'free')+
  guides(alpha = 'none') +
  labs(x = expression("R"^{2}))+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = 'grey50', linewidth = 0.2),
    strip.background = element_rect(colour = 'grey50', linewidth = 0.2)
  ); perm_out_full

ggsave('out/manuscript/supp/permanova.pdf', bg = 'white', 
       width = 2300, height = 1400, units = 'px', dpi = 280)
 
# # MONTHS TOGETHER -------
# 
# # 3 year models, bacteria 
# perms$Y3_16S <- adonis2(
#   betadiv$Bacteria$y3$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Time * Location, 
#   betadiv$Bacteria$y3$Meta, 
#   parallel = 6,
#   permutations = 999, method = "wunifrac", by = "terms"); perms$Y3_16S
# 
# # 3 year models, fungi
# perms$Y3_ITS <- adonis2(
#   betadiv$Fungi$y3$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z  + Year * Time * Location, 
#   betadiv$Fungi$y3$Meta, 
#   parallel = 6,
#   permutations = 999, method = "wunifrac", by = "terms"); perms$Y3_ITS
# 
# # 2-year models, bacteria
# perms$Y2_16S <- adonis2(
#   betadiv$Bacteria$y2$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Time * Location * Management, 
#   betadiv$Bacteria$y2$Meta, 
#   parallel = 6,
#   permutations = 999, method = "wunifrac", by = "terms"); perms$Y2_16S
# 
# # 2-year models, fungi
# perms$Y2_ITS <- adonis2(
#   betadiv$Fungi$y2$Dist$wuf ~ seq_depth_z + mean_temp + precip_72h_z + deg_h_z + Year * Time * Location * Management, 
#   betadiv$Fungi$y2$Meta, 
#   parallel = 6,
#   permutations = 999, method = "wunifrac", by = "terms"); perms$Y2_ITS
