pacman::p_load(phyloseq, tidyverse, mgx.tools, update = FALSE) # pak::pkg_install("jorondo1/mgx.tools?reinstall", upgrade = TRUE); 

ps.16S.ls <- read_rds('data/phyloseq/ps_filt_16S.rds')
ps.ITS.ls <- read_rds('data/phyloseq/ps_filt_ITS.rds')

# Diversity data, 
## species cluster aggolmeration
ps.list.out$Bacteria_sp_clust <-  tax_glom2(ps.list.out$Bacteria, taxrank = 'Species_cluster') 
ps.list.out$Fungi_sp_clust <-  tax_glom2(ps.list.out$Fungi, taxrank = 'Species_cluster') 

div.out <- list(
  Bacteria = mgx.tools::rarefy_diversity(
    ps.list.out$Bacteria, depth = 7000,
    n_iter=100, mc.cores = 8, vst = TRUE),
  
  Fungi =  mgx.tools::rarefy_diversity(
    ps.list.out$Fungi, depth = 3000, 
    n_iter=100, mc.cores = 8, vst = TRUE),
  
  Bacteria_sp_clust = mgx.tools::rarefy_diversity(
    ps.list.out$Bacteria_sp_clust, depth = 7000, 
    n_iter=100, mc.cores = 8, vst = TRUE),
  
  Fungi_sp_clust = mgx.tools::rarefy_diversity(
    ps.list.out$Fungi_sp_clust, depth = 3000, 
    n_iter=100, mc.cores = 8, vst = TRUE)
)

write_rds(div.out,
          file.path(analysis_shared_path,'R_data/reorganize/diversity_data.rds'), 
          compress = 'xz')
