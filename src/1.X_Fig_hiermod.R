# Alpha diversity figuress

pacman::p_load(tidyverse, purrr, patchwork)

## Setup -----------------------------------------------

pc_MDSTYCV_16S <- readRDS("out/hiermod/16S_7_tree_full_MDSTYCV/pc_MDSTYCV.rds")
pc_MDSTYCV_ITS <- readRDS("out/hiermod/ITS_7_tree_full_MDSTYCV/pc_MDSTYCV.rds")

source('src/utils/postcontrast_helpers.R')
source('src/hiermod/0_SETUP.R')

# Old plot version:
contrast_plot_panels(
  pc_MDSTYCV_16S$means, quant = c(0.001, 0.999), scales = 'free_y',
  group_pal = Management_palette,
  legend_title = "Posteriors (population medians)",
  ratio_stats = names(Fold_change_palette), ratio_pal = Fold_change_palette)


## Data wrangling --------------------

# Merge both datasets
plot_dat_raw <- rbind(
  pc_MDSTYCV_16S$means %>% mutate(Barcode = 'Bacteria'),
  pc_MDSTYCV_ITS$means %>% mutate(Barcode = 'Fungi')
) 

# Make contrasts Conventional - Organic
plot_dat <- plot_dat_raw %>% 
  mutate(
    value = case_when(
      group == 'Contrast' & str_detect(statistic, 'mean') ~ -value,
      TRUE~value
    ),
    Dataset = case_when(
      str_detect(statistic, 'fold difference') ~ 'Fold',
      TRUE ~ 'Diversity'
    )
  )

# and create column to split datasets (fold vs quantities)

# 89% Hpdi per group

summary_stats <- plot_dat %>% 
  group_by(Barcode, group, statistic) %>% 
  summarise(median = median(value),
            HPDI_lower = HPDI(value)[1],
            HPDI_upper = HPDI(value)[2],
            pd = max(mean(value > 0), mean(value<0)),
            .groups = 'drop')


# Unified legend colours (group + statistic variable)
combined_pal <- c(
  Management_palette[1:3], Fold_change_palette
)

# Management contrast plots

mean_plots <- plot_dat %>%
  filter(Dataset == 'Diversity',
         value <= quantile(value, 0.9999)) %>% 
  group_by(Dataset, Barcode) %>%
  group_split() %>%
  purrr::map(function(dat){
    
    dat %>% 
      ggplot(aes(x = value, fill = group, colour = group)) +
      geom_density(alpha = 0.5, linewidth = 0.1) +
      geom_vline(xintercept = 0, colour = 'grey60') +
      facet_grid(statistic~Barcode, scales = 'free') +
      labs(x = 'Effective number of ASVs')
  })

bact_means <- mean_plots[[1]] +
  theme(strip.text.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

fung_means <- mean_plots[[2]] +
  theme(axis.title.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

# Legend hidden on both mean panels 
mean_plot <- (bact_means + fung_means) &
  scale_fill_manual(values = combined_pal, limits = names(combined_pal), guide = "none") &
  scale_colour_manual(values = combined_pal, limits = names(combined_pal), guide = "none")


# Fold-diff plots
# scale_*_manual still draws legend keys for the Management-palette entries; 
# a legend key for a level absent from every plotted layer's data doesn't
# pick up the scale's assigned aesthetic).

fold_dummy <- tibble(group = names(Management_palette[1:3]), value = 1,
                     Barcode = plot_dat$Barcode[1])

# plot: 
fold_plots <- plot_dat %>%
  filter(Dataset == 'Fold') %>%
  mutate(group = as.character(statistic),
         strip_text = 'Fold changes') %>%
  group_by(Barcode) %>% 
  group_split() %>%
  purrr::map(function(dat){
    dat %>% 
      ggplot(bact_means$mapping) +
      geom_density(alpha = 0.5, linewidth = 0.1) +
      # dummy layer for legend:
      geom_area(data = fold_dummy, aes(x = value, y = 0, fill = group, colour = group),
                alpha = 0.5, linewidth = 0.1, inherit.aes = FALSE) +
      facet_grid(strip_text~., scales = 'free') +
      geom_vline(xintercept = 1, colour = "grey50") +
      labs(x = "Fold change in effective number of ASVs (Conventional/organic)", fill = NULL, colour = NULL) +
      # Force scale to show 0, 1, and integers only (no decimal)
      scale_x_continuous(
        breaks = function(x) {
          b <- scales::breaks_pretty(n = 4)(x)
          b <- unique(round(b))
          sort(unique(c(0, 1, b)))
        },
        limits = function(x) c(min(0, x[1]), x[2])  ) +
      theme(strip.text.x = element_blank())
  })

bact_fold <- fold_plots[[1]]  +
  theme(strip.text.y = element_blank())
fung_fold <- fold_plots[[2]] +  
  theme(axis.title.y = element_blank())

fold_plot <- (bact_fold + fung_fold)&
  scale_fill_manual(values = combined_pal, limits = names(combined_pal)) &
  scale_colour_manual(values = combined_pal, limits = names(combined_pal)) &
  theme(plot.margin = margin(20, 0, 0, 0))


# Patchwork:
adiv_plot <- patchwork::wrap_plots(
  list(mean_plot, fold_plot), 
  ncol = 1, heights = c(2,1)) +
  plot_layout(guides = 'collect') &
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = 'bottom',
    axis.title = element_text(size = 8)) & 
  labs(y = 'Density'); adiv_plot 


ggsave(plot = adiv_plot , filename = "out/manuscript/1_contrasts.pdf", bg = 'white', 
       width = 2400, height = 1800, units = 'px', dpi = 250)









