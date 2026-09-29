pacman::p_load(phyloseq, tidyverse, mgx.tools, update = FALSE) 
# pak::pkg_install("jorondo1/mgx.tools?reinstall", upgrade = TRUE); 

# Only Bacteria/Fungi -- the 2yr/3yr subsets in ps_objects_full.rds are
# derived from these two and don't need their own rarefaction.
ps.ls <- read_rds('data/ps_objects_full.rds')[c('Bacteria', 'Fungi')]

# Species-cluster agglomeration ---------------------------------------------
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
  Fungi = ps_fung,
  Bacteria_sp_clust = ps_bact_clust,
  Fungi_sp_clust = ps_fung_clust
)

write_rds(div.out,
          'data/diversity_data.rds',
          compress = 'xz')

# To patch in a metadata column without rerunning rarefaction, see
# src/archive/0.4.1_patch_diversity_dataset_column.R.
