# Checking correlation between read countst and diveristy indices-
# checking if that changes depending on rarefaction

pacman::p_load(phyloseq, tidyverse, ape, phangorn, btools)

ps.ls     <- read_rds('data/ps_objects_full.rds') # $Bacteria, $Fungi -- filtered, NOT rarefied
div_data  <- read_rds('data/diversity_data.rds')  # $Bacteria, $Fungi -- rarefied
ps.ls_sub <- list()
ps.ls_sub$Bacteria <- ps.ls$Bacteria
ps.ls_sub$Fungi <- ps.ls$Fungi

## Non rarefied diversity ------------------------------------------------------

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

nonrare <- map(ps.ls_sub, compute_nonrare_diversity) # $Bacteria, $Fungi

## Join against the rarefied indices ------------------------------------------

rare_cols <- c('Sample', 'Seq_depth', 'Richness', 'Shannon', 'Hill_1',
               'Simpson', 'Hill_2', 'Tail', 'Faith')

div_compare <- imap(nonrare, function(nr, kingdom){
  div_data[[kingdom]]$alpha %>%
    dplyr::select(any_of(rare_cols)) %>%
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
# its own panel 

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
    # overlay rarefied and non-rarefied
    scale_colour_manual(values = c(Rarefied = '#4E79A7', `Non-rarefied` = '#E15759')) +
    labs(x = 'Sequencing depth (raw, pre-rarefaction reads)', y = NULL, colour = NULL,
         title = 'Diversity index vs sequencing depth, rarefied vs non-rarefied') +
    theme_light() +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
  
})


# Spearman rank correlation with Seq_depth, per index/type, annotated onto
# its own panel 

#map(c('Fungi', 'Bacteria'), function(kingdom) {

cor_labels <- div_long %>%
  filter(index == 'Hill_1') %>% 
  group_by(Kingdom, type) %>%
  summarise(rho = cor(Seq_depth, value, method = 'spearman', use = 'complete.obs'), .groups = 'drop') %>%
  mutate(label = paste0(type, ': rho = ', round(rho, 2)),
         vjust = if_else(type == 'Rarefied', 1.5, 3))

p_rare_vs_nonrare <- div_long %>%
  filter(index == 'Hill_1') %>% 
  ggplot(aes(x = Seq_depth, y = value, colour = type)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_text(data = cor_labels, aes(x = Inf, y = Inf, label = label, vjust = vjust),
            hjust = 1.05, size = 3, show.legend = FALSE) +
  facet_grid(Kingdom~., scales = 'free') +
  # overlay rarefied and non-rarefied
  scale_colour_manual(values = c(Rarefied = '#4E79A7', `Non-rarefied` = '#E15759')) +
  labs(x = 'Samlpe read count', y = 'Effective number of ASVs (Hill 1)', colour = NULL) +
  theme(legend.position = c(0.8,0.8),
        legend.background = element_rect(
          linewidth = 0.2,
          colour = 'black'
        )); p_rare_vs_nonrare

ggsave(paste0('out/summaries/seqdepth_cor.pdf'),
       p_rare_vs_nonrare, width = 2000, height = 1300,
       units = 'px', dpi = 280, bg = 'white')
# })


