library(phyloseq)
library(dplyr)
library(ggplot2)
library(patchwork)

###Seperate Data into the 2 datasets

ps.ls <- ps_objects_full

cultivars_3yr <- c("Cortland", "Liberty", "Paulared")
cultivars_2yr <- c("Honeycrisp", "Spartan")

build_subset_ps <- function(ps, keep_cultivars) {
  
  # Get metadata
  df <- data.frame(
    sample_data(ps),
    stringsAsFactors = FALSE
  )
  
  # Identify samples belonging to desired cultivars
  keep_samples <- rownames(df)[
    df$Cultivar %in% keep_cultivars
  ]
  
  # Keep only those samples
  ps_sub <- prune_samples(
    keep_samples,
    ps
  )
  
  # Remove ASVs that have zero reads after subsetting
  ps_sub <- prune_taxa(
    taxa_sums(ps_sub) > 0,
    ps_sub
  )
  
  return(ps_sub)
}

Years3 <- lapply(ps.ls, build_subset_ps, keep_cultivars = cultivars_3yr)
Years2 <- lapply(ps.ls, build_subset_ps, keep_cultivars = cultivars_2yr)

names(Years3) <- names(ps.ls)
names(Years2) <- names(ps.ls)

sapply(ps.ls, nsamples)
sapply(Years3, nsamples)
sapply(Years2, nsamples)

saveRDS(Years3, file = "Years3.rds")
saveRDS(Years2, file = "Years2.rds")