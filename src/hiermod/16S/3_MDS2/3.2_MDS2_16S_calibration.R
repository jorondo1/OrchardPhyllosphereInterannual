# 16S, SHIFTED (Hill_1 - 1): Management x Season

# Location REMOVED because it is too confounded with Management
# Willdo sensitivity test later on with an unconfounded (though less powerful)
# subset.

# sigma stays Management-only (no cell-level split yet). Parameter
# recovery, prior-predictive check, and SBC via the SBC package -- same
# template as MDL's own rewritten 2.2 script.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MDS2_model.R') # model_MDS2_16S, means_MDS2(), sim_div_MDS2(), simulate_from_priors_MDS2()

model <- model_MDS2_16S

hiermod_out_dir <- "out/hiermod/16S_3_lognormal_MDS2_shifted"

## Parameter recovery -----------------------------------------------------------
# Same baseline/gap values as MDL's and the old MDLS2 family's own
# calibration, so results stay comparable.

may_conv <- 180
may_org  <- 120
july_conv_shift <- 0.35
july_org_shift  <- 0.2

true_sigma <- cv_to_sigma(c(0.5, 0.8)) # conv, org

set.seed(20260911)

dat_sim <- sim_div_MDS2(
  N_samples = 240,
  loga = log(c(may_conv, may_org)),
  s_conv = july_conv_shift,
  gap_shift = july_org_shift,
  sigma = true_sigma,
  year_offset = c(0, 0.3, -0.2),
  sigma_tr = 0.3,
  shift = 1
); head(dat_sim)

hist(dat_sim$Dv_shifted, breaks = 100)

fit_sim <- ulam(
  model,
  data = as.list(dat_sim),
  chains = 6, cores = 6, iter = 5000,
  control = list(adapt_delta = 0.99))
precis(fit_sim, depth = 2)

post_sim <- extract.samples(fit_sim)

### Fixed effect + sigma recovery ---------

(param_recovery <- check_recovery(
  true = list(
    loga1 = log(may_conv), loga2 = log(may_org),
    sigma1 = true_sigma[1], sigma2 = true_sigma[2]),
  post_draws = list(
    loga1 = post_sim$loga[,1], loga2 = post_sim$loga[,2],
    sigma1 = post_sim$sigma[,1], sigma2 = post_sim$sigma[,2])))
### Contrast recovery -------------

cr <- contrast_recovery(
  fit_sim, means_MDS2, may_conv, may_org, july_conv_shift, july_org_shift, shift = 1)

p_sim_contrast <- contrast_plot_panels(
  cr$estimands, quant = c(0.01, 0.99), scales = 'free_y',
  group_pal = Management_palette,
  true_vals = cr$true_estimands); p_sim_contrast

save_report("sim_summary", model_id, recovery = param_recovery, fit_sim, cr$estimands, model, model_name = "The Locationless")
save_gg("sim_contrast_density", model_id, p_sim_contrast)

save_pdf("sim_trankplot", model_id, 
         function() trankplot(fit_sim, max_rows = 30, n_cols = 5),
         width = 30, height = 50)

## Prior predictive check -------------------------------------------------------

n_prior <- 1000
extracted_prior <- extract.prior(fit_sim, n = n_prior)

prior_pred <- map_dfr(seq_len(n_prior), function(i){
  simulate_from_priors_MDS2(draw_true(extracted_prior, i), shift = 1)
}, .id = "draw")

summary(prior_pred$Dv_shifted) # judge on median/IQR, not mean/SD

p_prior_pc <- prior_predictive_spaghetti(
  prior_pred, value_col = "Dv_shifted", upper_q = 0.99, model = model,
  title = "Prior predictive check", observed = dat_sim$Dv_shifted); p_prior_pc

save_gg("sim_prior_PC", model_id, p_prior_pc)

## loga x sigma[Mg]/sigma_tr funnel check ------------------------------------
# SBC on the real 100-replicate run found loga[]'s miscalibration is just as
# severe here (no Location at all) as it was in MDLS2v (Location present) --
# ruling out Location's low cardinality as the driver. sigma[Mg] is now the
# only remaining scale parameter loga could be entangled with in the same
# likelihood term (Dv ~ dlnorm(mu, sigma[Mg])); sigma_tr checked too since it
# showed a secondary correlation with the treedepth blowups.

p_funnel <- function(){
  par(mfrow = c(2,2))
  plot(post_sim$sigma[,1], post_sim$loga[,1],
       xlab = "sigma[1] (Conventional)", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma[,2], post_sim$loga[,2],
       xlab = "sigma[2] (Organic)", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma_tr, post_sim$loga[,1],
       xlab = "sigma_tr", ylab = "loga[1] (Conventional)", pch = 16, col = scales::alpha("black", 0.15))
  plot(post_sim$sigma_tr, post_sim$loga[,2],
       xlab = "sigma_tr", ylab = "loga[2] (Organic)", pch = 16, col = scales::alpha("black", 0.15))
  par(mfrow = c(1,1))
}
save_pdf("loga_sigma_funnel", model_id, p_funnel)

## Simulation-based calibration (SBC), via the SBC package -----------------------
# Same template as MDL's/MDLS2v's own calibration. yr[Yr]/tr[Tr] stay out
# of `variables` for the same reason as those models: the simulator draws
# fresh per-level offsets internally, so a prior draw of the array isn't
# what generated the data.
sbc_gen_MDS2 <- make_sbc_generator(
  fit = fit_sim, simulate_fn = simulate_from_priors_MDS2,
  keep = c("loga", "s_conv", "gap_shift", "sigma", "sigma_tr"),
  gen_cols = c("Dv", "Mg", "Mo", "Yr", "Tr"),
  extra_globals = "sim_div_MDS2", shift = 1)

n_sbc  <- 100
n_iter <- 10000

sbc_MDS2 <- run_sbc_pipeline(
  generator = sbc_gen_MDS2$generator, globals = sbc_gen_MDS2$globals,
  n_sbc = n_sbc, model = model, model_id = model_id, n_iter = n_iter,
  hiermod_out_dir = hiermod_out_dir, dquants = dq_MDS2,
  control = list(adapt_delta = 0.99))

plot_sbc_diagnostics(sbc_MDS2, model_id, n_sbc)

sbc_MDS2$stats |>
  dplyr::filter(variable %in% c("loga[1]", "loga[2]")) |>
  dplyr::group_by(variable) |>
  dplyr::summarise(mean_rank_frac = mean(rank / max_rank), median_rank_frac = median(rank / max_rank))

save_sbc_health_report(model_id, sbc_MDS2, n_sbc, n_iter,
                        variables = c("loga[1]", "loga[2]", "sigma[1]", "sigma[2]", "sigma_tr",
                                      "may_gap", "july_gap", "seasonal_change"),
                        hiermod_out_dir = hiermod_out_dir)
