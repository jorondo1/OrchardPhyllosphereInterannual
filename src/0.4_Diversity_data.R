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
