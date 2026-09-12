# MODEL 5 (MDLS2), 16S: prior sensitivity check. Refits with looser sigma
# priors (dexp(1.5) instead of the calibrated dexp(2.5)) and compares 
# posteriors against real fit. If they barely move, the calibrated
# prior isn't the thing driving the result.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDLS2_model.R')

hiermod_out_dir <- "out/hiermod/16S_5_lognormal_MDLS2_shifted"
fit_MDLS2  <- readRDS(file.path(hiermod_out_dir, "fit_MDLS2_shifted.rds"))
dat        <- readRDS(file.path(hiermod_out_dir, "dat_MDLS2_shifted.rds"))
model_ppc1 <- readRDS(file.path(hiermod_out_dir, "model_ppc1.rds"))

model_loose <- model_ppc1
model_loose$pr_sigma     <- quote(sigma[cell] ~ dexp(1.5))
model_loose$pr_sigma_loc <- quote(sigma_loc   ~ dexp(1.5))
model_loose$pr_sigma_tr  <- quote(sigma_tr    ~ dexp(1.5))

fit_loose <- ulam(
  model_loose,
  data = dat,
  chains = 6, cores = 6, iter = 20000,
  control = list(adapt_delta = 0.99)
)
save_fit("fit", "MDLS2_shifted_loose", fit_loose)

## Compare posteriors -----------------------------------------------

pf_tight <- post_full(fit_MDLS2, means_MDLS2, shift = 1)
pf_loose <- post_full(fit_loose, means_MDLS2, shift = 1)

summarize_fit <- function(pf, label){
  tibble(
    fit           = label,
    sigma1        = median(pf$sigma$sigma_1),
    sigma2        = median(pf$sigma$sigma_2),
    sigma3        = median(pf$sigma$sigma_3),
    sigma4        = median(pf$sigma$sigma_4),
    sigma_loc     = median(pf$sigma_loc$sigma_loc),
    sigma_tr      = median(pf$sigma_tr$sigma_tr),
    may_contrast  = median(pf$mean$mean_3 - pf$mean$mean_1),
    july_contrast = median(pf$mean$mean_4 - pf$mean$mean_2)
  )
}

sensitivity_comparison <- bind_rows(
  summarize_fit(pf_tight, "tight (dexp(2.5) for all sigmas)"),
  summarize_fit(pf_loose, "loose (dexp(1.5) for all sigmas)")
)
print(sensitivity_comparison)
# parameters barely move except sigma_loc from .274 to .300, which seems minor
write.csv(sensitivity_comparison, file.path(hiermod_out_dir, "sensitivity_comparison.csv"), row.names = FALSE)

