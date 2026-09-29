# PERMANOVA/betadisper tests (weighted UniFrac) behind the main+supp PCoA
# figures, both barcodes: dispersion homogeneity per grouping variable, and
# the two adonis2() formulas (with/without environmental covariates), for
# 2y/3y and their May/July splits.

pacman::p_load(tidyverse, vegan, update = FALSE)

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
  betadiv$Bacteria$y3_may$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Location, 
  betadiv$Bacteria$y3_may$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["3Y - July 16S"]] <- adonis2(
  betadiv$Bacteria$y3_july$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Location, 
  betadiv$Bacteria$y3_july$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["2Y - May 16S"]] <- adonis2(
  betadiv$Bacteria$y2_may$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Location * Management, 
  betadiv$Bacteria$y2_may$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["2Y - July 16S"]] <- adonis2(
  betadiv$Bacteria$y2_july$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Location * Management, 
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
  betadiv$Fungi$y3_may$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Location, 
  betadiv$Fungi$y3_may$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["3Y - July ITS"]] <- adonis2(
  betadiv$Fungi$y3_july$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp  + Year * Location, 
  betadiv$Fungi$y3_july$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["2Y - May ITS"]] <- adonis2(
  betadiv$Fungi$y2_may$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Location * Management, 
  betadiv$Fungi$y2_may$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")

perms[["2Y - July ITS"]] <- adonis2(
  betadiv$Fungi$y2_july$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Location * Management, 
  betadiv$Fungi$y2_july$Meta, 
  parallel = 6,
  permutations = 999, method = "wunifrac", by = "terms")


# Export: one sheet per perms element ------
save_stat_xlsx("out/manuscript/betadiv_permanova.xlsx", perms)

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
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

ggsave('out/manuscript/supp/permanova.pdf', bg = 'white', 
       width = 2300, height = 1700, units = 'px', dpi = 200)

# 
# 
# # MONTHS TOGETHER -------
# 
# # 3 year models, bacteria 
# perms$Y3_16S <- adonis2(
#   betadiv$Bacteria$y3$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Time * Location, 
#   betadiv$Bacteria$y3$Meta, 
#   parallel = 6,
#   permutations = 999, method = "wunifrac", by = "terms"); perms$Y3_16S
# 
# # 3 year models, fungi
# perms$Y3_ITS <- adonis2(
#   betadiv$Fungi$y3$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp  + Year * Time * Location, 
#   betadiv$Fungi$y3$Meta, 
#   parallel = 6,
#   permutations = 999, method = "wunifrac", by = "terms"); perms$Y3_ITS
# 
# # 2-year models, bacteria
# perms$Y2_16S <- adonis2(
#   betadiv$Bacteria$y2$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Time * Location * Management, 
#   betadiv$Bacteria$y2$Meta, 
#   parallel = 6,
#   permutations = 999, method = "wunifrac", by = "terms"); perms$Y2_16S
# 
# # 2-year models, fungi
# perms$Y2_ITS <- adonis2(
#   betadiv$Fungi$y2$Dist$wuf ~ seq_depth_z + deg_h_z + precip_72h_z + mean_temp + Year * Time * Location * Management, 
#   betadiv$Fungi$y2$Meta, 
#   parallel = 6,
#   permutations = 999, method = "wunifrac", by = "terms"); perms$Y2_ITS
