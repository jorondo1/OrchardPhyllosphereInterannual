# Exploratory, to be archived
# - 1: env. variable vs PCoA axis correlations + residual checks
# - 2: PERMANOVA blocked by tree

pacman::p_load(tidyverse, vegan, patchwork, update = FALSE)  # patchwork: scatter + qq below

source('src/utils/beta_div_saver.R') # save_stat_kable(), mm_to_in()

out_dir <- "out/exploration/betadiv/envfit_corr"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(230726)

betadiv <- readRDS('data/diversity_subsets.rds')

## PCoA coordinates + metadata, every dataset x metric -------------------------------
pcoa_data <- expand_grid(
  kingdom = c("Bacteria", "Fungi"), dataset = c("2Y", "3Y"),
  metric = c("wuf", "bray")) %>%
  mutate(
    label  = paste0(dataset, substr(kingdom, 1, 1)),  # "2YB", "3YB", "2YF", "3YF"
    dist   = pmap(list(kingdom, dataset, metric), ~ betadiv[[..1]][[tolower(paste0("y", substr(..2, 1, 1)))]]$Dist[[..3]]),
    meta   = map2(kingdom, dataset, ~ betadiv[[.x]][[tolower(paste0("y", substr(.y, 1, 1)))]]$Meta),
    pcoa_df = map2(dist, meta, function(d, m) {
      pcoa <- capscale(d ~ 1)
      scores(pcoa, display = "sites") %>% as.data.frame() %>%
        rownames_to_column("Sample") %>%
        select(Sample, PCo1 = MDS1, PCo2 = MDS2) %>%
        left_join(m, by = "Sample")
    }))

## Correlation diagnostics for one env. variable x one PCoA axis, z-scored --------------
# TODO: uses raw mean_temp/precip_72h/deg_h, not the centered *_z versions
alpha <- 0.05  # significance threshold, lm slope + Shapiro-Wilk test

# Correlation of a z-scored variable with one PCoA axis + residual normality checks
corr_diagnostics <- function(df, var, axis) {
  axis_value <- df[[axis]]
  z <- as.numeric(scale(df[[var]]))
  model <- lm(z ~ axis_value)

  r <- cor(z, axis_value, use = "complete.obs")
  p <- summary(model)$coefficients[2, "Pr(>|t|)"]
  shapiro_p <- stats::shapiro.test(residuals(model))$p.value

  scatter <- ggplot(data.frame(axis_value, z), aes(x = axis_value, y = z)) +
    geom_point(colour = "grey40", alpha = 0.6) +
    geom_smooth(method = "lm", colour = "black", fill = "grey70", linewidth = 0.7) +
    theme_bw() +
    labs(title = "Standardized relationship", x = axis, y = paste0(var, " (z-score)"))

  qq <- ggplot(data.frame(resid = residuals(model)), aes(sample = resid)) +
    stat_qq(colour = "grey40") +
    stat_qq_line(colour = "black") +
    theme_bw() +
    labs(title = "Residual Q-Q plot", x = "Theoretical quantiles", y = "Residual quantiles")

  list(r = r, p = p, shapiro_p = shapiro_p,
       significant = p < alpha, normal = shapiro_p > alpha,
       diagnostic_plot = scatter + qq)
}

env_vars <- c("mean_temp", "precip_72h", "deg_h")
axes     <- c("PCo1", "PCo2")

## Correlation grid: variable x axis, for every dataset x metric ------------------------
corr_results <- pcoa_data %>%
  mutate(grid = map(pcoa_df, ~ tidyr::expand_grid(variable = env_vars, axis = axes))) %>%
  tidyr::unnest(grid) %>%
  mutate(diagnostics = pmap(list(pcoa_df, variable, axis), corr_diagnostics))

## One summary table: r/p/shapiro_p/significant/normal, per dataset x variable x axis x metric --
corr_diag_summary <- corr_results %>%
  mutate(
    r = map_dbl(diagnostics, "r"),
    p = map_dbl(diagnostics, "p"),
    shapiro_p = map_dbl(diagnostics, "shapiro_p"),
    significant = map_lgl(diagnostics, "significant"),
    normal = map_lgl(diagnostics, "normal")
  ) %>%
  select(metric, label, variable, axis, r, p, shapiro_p, significant, normal)

save_stat_kable(
  file.path(out_dir, "CorrDiagnostics_summary.html"), corr_diag_summary, rownames_to = NULL,
  caption = paste0("Environmental-variable correlation diagnostics: r/p from lm(z-scored variable ~ PCoA axis); ",
                    "shapiro_p from a Shapiro-Wilk test on the model residuals (significant = lm slope p < alpha; ",
                    "normal = Shapiro-Wilk p > alpha).")
)

## One multi-page PDF per dataset x metric: scatter+QQ plot pair, one pair per page --
purrr::walk(unique(corr_results$metric), function(met) {
  purrr::walk(unique(corr_results$label), function(lab) {
    plots <- corr_results %>%
      filter(metric == met, label == lab) %>%
      pull(diagnostics) %>% map("diagnostic_plot")
    path <- file.path(out_dir, paste0("CorrDiagnostics_", lab, "_", met, ".pdf"))
    grDevices::pdf(path, width = mm_to_in(170), height = mm_to_in(90))
    purrr::walk(plots, print)
    grDevices::dev.off()
  })
})

## Section 2: PERMANOVA, blocked by tree -----------------------------------------------
# Code nested within Management. One HTML table per dataset x metric.
pcoa_data %>%
  mutate(adonis = map2(dist, pcoa_df, ~ adonis2(
    .x ~ Year * Time * (Management / Code), data = .y,
    permutations = how(nperm = 999, blocks = .y$Tree_id)))) %>%
  { walk(seq_len(nrow(.)), function(i) save_stat_kable(
      file.path(out_dir, paste0("Permanova_", .$label[i], "_", .$metric[i], ".html")), .$adonis[[i]],
      caption = paste0(.$label[i], " (", .$metric[i], "): PERMANOVA (adonis2), blocked by Tree_id"))) }
