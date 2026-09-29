# Beta-diversity data prep: one nested list holding metadata + distance
# matrices for every subset used by the fig/stats/envfit scripts.
#   betadiv$<Barcode>$<subset>$Meta          -- metadata tibble
#   betadiv$<Barcode>$<subset>$Dist$<metric> -- dist object (bray/wuf)

pacman::p_load(tidyverse, update = FALSE)

div_obj <- readRDS("data/diversity_data.rds")

# Slice a dist object down to one sample set ---------------------------------
subset_dist <- function(dist_obj, samples) {
  mat <- as.matrix(dist_obj)
  keep <- rownames(mat) %in% samples
  as.dist(mat[keep, keep])
}

# Meta + Dist for one sample set ---------------------------------------------
make_entry <- function(meta, beta) {
  list(
    Meta = meta,
    Dist = list(
      bray = subset_dist(beta$bray, meta$Sample),
      wuf  = subset_dist(beta$unifrac_w, meta$Sample)
    )
  )
}

# All 7 subsets for one barcode ----------------------------------------------
build_subsets <- function(alpha, beta) {
  meta_all <- alpha %>% mutate(Time = factor(Time, levels = c("May", "July")))
  meta_y2  <- meta_all %>% filter(Dataset == "2-year")
  meta_y3  <- meta_all %>% filter(Dataset == "3-year")

  list(
    all     = make_entry(meta_all, beta),
    y2      = make_entry(meta_y2, beta),
    y2_may  = make_entry(filter(meta_y2, Time == "May"), beta),
    y2_july = make_entry(filter(meta_y2, Time == "July"), beta),
    y3      = make_entry(meta_y3, beta),
    y3_may  = make_entry(filter(meta_y3, Time == "May"), beta),
    y3_july = make_entry(filter(meta_y3, Time == "July"), beta)
  )
}

betadiv <- list(
  Bacteria = build_subsets(div_obj$Bacteria$alpha, div_obj$Bacteria$beta),
  Fungi    = build_subsets(div_obj$Fungi$alpha, div_obj$Fungi$beta)
)

saveRDS(betadiv, 'data/diversity_subsets.rds', compress = 'xz')

