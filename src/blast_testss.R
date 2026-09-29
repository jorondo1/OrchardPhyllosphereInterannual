pacman::p_load(tidyverse, Biostrings, phyloseq)


# load
ps.ls <- read_rds('data/ps_objects_full.rds')

# extract tax table as df
asvs<- ps.ls$Fungi %>% 
  tax_table() %>% 
  as.data.frame() %>% 
  rownames_to_column('ASV')


# helper function for a single species cluster
blast_fun <- function(asvs, name) {
  

  dna <- asvs %>% filter(Species_cluster==name)
  out <- dna %>% pull(ASV) %>% Biostrings::DNAStringSet()
  names(out) <- dna %>% pull(Species_cluster)

  out  
}

# single use 
blast_fun(asvs, 'Helotiales_sp_clust_17')


# list to loop; 
sp_clust_to_check <- list(
 # 'Helotiales_sp_clust_17',
#  'Setomelanomma_sp_clust_23',
  'Ascomycota_sp_clust_5',
  'NA_sp_clust_4',
  'Pleosporales_gen_Incertae_sedis_sp_clust_7',
  'Ascomycota_sp_clust_14',
  'NA_sp_clust_15'
  )

out <- map(sp_clust_to_check, ~ blast_fun(asvs, .x)) %>% 
  DNAStringSetList() %>% 
  unlist()

Biostrings::writeXStringSet(out, paste0('out/exploration/asvs_to_blast.fasta'))




