# Checking rarefaction 
pacman::p_load(phyloseq, tidyverse, ape, phangorn, btools)

ps.ls    <- read_rds('data/ps_objects_full.rds') # $Bacteria, $Fungi -- filtered, NOT rarefied
div.ITS  <- read_rds('data/diversity_data.rds')

## Non-rarefied diversity ------------------------------------------------------
# Same formulas as mgx.tools::rarefy_diversity()/compute_diversities.R,
# applied once directly to the raw count table instead of averaging over
# repeated rarefaction draws. Richness/Shannon/Simpson (and the two Hill
# numbers built from them) use relative abundances, matching the rarefied
# calculation; Tail is reproduced from raw counts (not proportions) because
# that's how it's defined upstream too, not a normalization we're skipping.

compute_nonrare_diversity <- function(ps){
  if (!ape::is.rooted(phy_tree(ps))) {
    phy_tree(ps) <- phangorn::midpoint(phy_tree(ps))
  }
  
  counts <- as(otu_table(ps), 'matrix')
  if (taxa_are_rows(ps)) counts <- t(counts) # samples x taxa
  
  row_sums <- rowSums(counts)
  richness <- rowSums(counts > 0)
  
  p       <- counts / row_sums
  logp    <- ifelse(p > 0, log(p), 0)
  shannon <- -rowSums(p * logp)
  simpson <- rowSums(p^2)
  
  n_taxa          <- ncol(counts)
  tail_multiplier <- (seq_len(n_taxa) - 1)^2
  sorted_mat      <- apply(counts, 1, sort, decreasing = TRUE) # n_taxa x n_samples
  tail_vals       <- sqrt(colSums(sorted_mat * tail_multiplier))
  
  faith_res <- suppressWarnings(suppressMessages(btools::estimate_pd(ps)))
  
  tibble(
    Sample           = rownames(counts),
    Richness_nonrare = richness,
    Shannon_nonrare  = shannon,
    Hill_1_nonrare   = exp(shannon),
    Simpson_nonrare  = simpson,
    Hill_2_nonrare   = 1 / simpson,
    Tail_nonrare     = tail_vals
  ) %>%
    left_join(
      faith_res %>% tibble::rownames_to_column('Sample') %>%
        transmute(Sample, Faith_nonrare = PD),
      by = 'Sample'
    )
}

nonrare <- map(ps.ls, compute_nonrare_diversity) # $Bacteria, $Fungi

## Join against the rarefied indices ------------------------------------------

rare_cols <- c('Sample', 'Seq_depth', 'Richness', 'Shannon', 'Hill_1',
               'Simpson', 'Hill_2', 'Tail', 'Faith')

div_compare <- imap(nonrare, function(nr, kingdom){
  div.ITS[[kingdom]]$alpha %>%
    select(any_of(rare_cols)) %>%
    left_join(nr, by = 'Sample')
}) %>% list_rbind(names_to = 'Kingdom')

## Long format for plotting ----------------------------------------------------

indices <- c('Richness', 'Shannon', 'Hill_1', 'Simpson', 'Hill_2', 'Tail', 'Faith')

div_long <- indices %>%
  map(function(idx){
    div_compare %>%
      transmute(Kingdom, Seq_depth,
                index = idx,
                Rarefied       = .data[[idx]],
                `Non-rarefied` = .data[[paste0(idx, '_nonrare')]]) %>%
      pivot_longer(c(Rarefied, `Non-rarefied`), names_to = 'type', values_to = 'value')
  }) %>%
  list_rbind() %>%
  mutate(index = factor(index, levels = indices))

## Plot -------------------------------------------------------------------------

# Spearman rank correlation with Seq_depth, per index/type, annotated onto
# its own panel -- text colour inherits the same `type` mapping as the
# points (aes(colour = type) is set once, at the top level), so no separate
# colour scale to keep in sync.

map(c('Fungi', 'Bacteria'), function(kingdom) {
  cor_labels <- div_long %>%
    filter(Kingdom == kingdom) %>%
    group_by(index, type) %>%
    summarise(rho = cor(Seq_depth, value, method = 'spearman', use = 'complete.obs'), .groups = 'drop') %>%
    mutate(label = paste0(type, ': rho = ', round(rho, 2)),
           vjust = if_else(type == 'Rarefied', 1.5, 3))
  p_rare_vs_nonrare <- div_long %>%
    filter(Kingdom == kingdom) %>%
    ggplot(aes(x = Seq_depth, y = value, colour = type)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_text(data = cor_labels, aes(x = Inf, y = Inf, label = label, vjust = vjust),
              hjust = 1.05, size = 3, show.legend = FALSE) +
    facet_grid(index ~ Kingdom, scales = 'free') +
    scale_colour_manual(values = c(Rarefied = '#4E79A7', `Non-rarefied` = '#E15759')) +
    labs(x = 'Sequencing depth (raw, pre-rarefaction reads)', y = NULL, colour = NULL,
         title = 'Diversity index vs sequencing depth, rarefied vs non-rarefied') +
    theme_light() +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
  
  ggsave(paste0('out/summaries/seqdepth_cor_',kingdom,'.pdf'),
         p_rare_vs_nonrare, width = 2000, height = 3000,
         units = 'px', dpi = 220, bg = 'white')
})
