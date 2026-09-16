# check_seq_depth_confound.R -- standalone diagnostic: does sequencing
# depth act as a confound between Management (and/or Month) and Hill_1
# diversity? Supports the "Control variables" TODO section
# (src/hiermod/16S/TODO.md) written after the hiermod MDSYC/MDSYCV models
# found the raw Management gap roughly halves once seq_depth_z (and
# weather) are controlled for -- that work was 16S-only, this check covers
# both Bacteria and Fungi (faceted by Barcode), not just 16S.
#
# Seq_depth here is the DADA2 pipeline's post-filtering (pre-rarefaction)
# read count, it's a proxy for true sequencing depth, not a direct
# measurement. 

pacman::p_load(tidyverse, rstatix, patchwork)

management_pal <- c(Conventional = "#F28E2B", Organic = "#499894")
month_pal      <- c(May = "#4E79A7", July = "#E15759")

div.ls <- read_rds('data/diversity_data.rds')
div <- rbind(
  div.ls$Bacteria$alpha %>% mutate(Barcode='Bacteria'),
  div.ls$Fungi$alpha %>% mutate(Barcode='Fungi')
) %>% group_by(Barcode)

## Panel 1: overall Hill_1 vs Seq_depth relationship --------------------------

ct <- div %>% 
  cor_test(vars = c('Hill_1', 'Seq_depth')) %>%
  mutate(
    label = sprintf(
      "r = %.2f, p %s",
      cor,
      ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
    )
  )

p1 <- ggplot(div, aes(x = Seq_depth, y = Hill_1)) +
  geom_point(alpha = 0.4, size = 1.5, colour = "grey30") +
  geom_smooth(method = "lm", colour = "black", se = TRUE) +
  facet_wrap(~Barcode, scales = 'free') +
  geom_text(
    data = ct,                       
    aes(x = Inf, y = Inf, label = label),  
    hjust = 1.1, vjust = 1.5, size = 4,
    inherit.aes = FALSE             
  ) +
  labs(subtitle = "Read count predicts Hill numbers",
       x = "Read count", y = "Hill number of order 1"); p1

## Panel 2: same relationship, by Management -----------------------------------
# Shared slope, different typical depth -- both pieces of the confound at
# once (deeper sequencing in Organic, per the t-test in the TODO entry).

# detailed = TRUE -- t_test() otherwise drops estimate/conf.low/conf.high,
# which the label needs. group1/group2 (here "Conventional"/"Organic", the
# levels' own alphabetical order) come straight from rstatix, so the label
# text can't drift out of sync with what `estimate` actually is.
mg_test <- div %>%
  t_test(Seq_depth ~ Management, detailed = TRUE) %>%
  mutate(
    label = sprintf(
      "Difference in mean read count\n%s - %s:\n%.0f [%.0f, %.0f],\nt = %.2f, p %s",
      group1, group2, estimate, conf.low, conf.high, statistic,
      ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
    )
  )

p2 <- ggplot(div, aes(x = Seq_depth, y = Hill_1, colour = Management)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_colour_manual(values = management_pal) +
  facet_wrap(~Barcode, scales = 'free') +
  geom_text(
    data = mg_test,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.02, vjust = 1.5, size = 4, colour = "grey20",
    inherit.aes = FALSE
  ) +
  labs(subtitle = "Relationship by management",
       x = "Read count", y = "Hill number of order 1", colour = NULL); p2

## Panel 3: same relationship, by Month ----------------------------------------
# "Are we seeing something alike between months?" -- yes, and for weather
# covariates even more so (see TODO entry); this panel checks Seq_depth
# itself.

mo_test <- div %>%
  t_test(Seq_depth ~ Time, detailed = TRUE) %>%
  mutate(
    label = sprintf(
      "Difference in mean read count\n%s - %s:\n%.0f [%.0f, %.0f],\nt = %.2f, p %s",
      group1, group2, estimate, conf.low, conf.high, statistic,
      ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
    )
  )

p3 <- ggplot(div, aes(x = Seq_depth, y = Hill_1, colour = Time)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_colour_manual(values = month_pal) +
  facet_wrap(~Barcode, scales = 'free') +
  geom_text(
    data = mo_test,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.02, vjust = 1.5, size = 4, colour = "grey20",
    inherit.aes = FALSE
  ) +
  labs(subtitle = "Relationshp by month",
       x = "Read count", y = "Hill number of order 1", colour = NULL); p3

p_confound <- (p1 / p2 / p3) & theme_light()
p_confound

dir.create("out/summaries", recursive = TRUE, showWarnings = FALSE)
ggsave("out/summaries/seqdepth_management_confound.pdf", p_confound,
       width = 2500, height = 3200, units = "px", dpi = 220, bg = "white")
