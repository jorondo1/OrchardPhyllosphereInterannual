pacman::p_load(phyloseq, tidyverse, mgx.tools, update = FALSE) # pak::pkg_install("jorondo1/mgx.tools?reinstall", upgrade = TRUE); 

ps.ls <- read_rds('data/ps_objects_full.rds')

# Diversity data, 
## species cluster aggolmeration
ps_clust.ls <- map(ps.ls, tax_glom2, taxrank = 'Species_cluster') 

ps_bact <- mgx.tools::rarefy_diversity(
  ps.ls$Bacteria, depth = 7000,
  n_iter=100, mc.cores = 7, vst = TRUE)

ps_fung <- mgx.tools::rarefy_diversity(
  ps.ls$Fungi, depth = 3000, 
  n_iter=100, mc.cores = 7, vst = TRUE)

ps_bact_clust <- mgx.tools::rarefy_diversity(
  ps_clust.ls$Bacteria, depth = 7000, 
  n_iter=100, mc.cores = 7, vst = TRUE)

ps_fung_clust <- mgx.tools::rarefy_diversity(
  ps_clust.ls$Fungi, depth = 3000, 
  n_iter=100, mc.cores = 7, vst = TRUE)

div.out <- list(
  Bacteria = ps_bact,
  
  Fungi =  ps_fung,
  
  Bacteria_sp_clust = ps_bact_clust,
  
  Fungi_sp_clust = ps_fung_clust
)

write_rds(div.out,
          'data/diversity_data.rds', 
          compress = 'xz')

# update metadata if need be
# div.out <- read_rds('data/diversity_data.rds')
# bact_samdat <- mgx.tools::samdat_as_tibble(ps.ls.out$Bacteria)
# fung_samdat <- mgx.tools::samdat_as_tibble(ps.ls.out$Fungi)
# 
# div.out$Fungi_sp_clust$alpha[, names(fung_samdat)] <- fung_samdat[, names(fung_samdat)]
# div.out$Fungi$alpha[, names(fung_samdat)]          <- fung_samdat[, names(fung_samdat)]
# 
# div.out$Bacteria_sp_clust$alpha[, names(bact_samdat)] <- bact_samdat[, names(bact_samdat)]
# div.out$Bacteria$alpha[, names(bact_samdat)]          <- bact_samdat[, names(bact_samdat)]
# 
# div.out$Bacteria_sp_clust$alpha %<>% 
#   left_join(bact_samdat %>% select(Sample, Seq_depth), 
#           by = 'Sample')
# 
# div.out$Bacteria$alpha %<>% 
#   left_join(bact_samdat %>% select(Sample, Seq_depth), 
#             by = 'Sample')
# 
# div.out$Fungi_sp_clust$alpha %<>% 
#   left_join(fung_samdat %>% select(Sample, Seq_depth), 
#             by = 'Sample')

write_rds(div.out,
          'data/diversity_data.rds', 
          compress = 'xz')
