# MODEL 1 (MD) + MODEL 2 (MDv), 16S: posterior contrast, run against the saved fits.

hiermod_marker <- "16S"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MD_16S/model_MDv_16S, means_MD/means_MDv

model_MD  <- model_MD_16S
model_MDv <- model_MDv_16S

# Updated priors (see MDv calibration)
model_MDv$prior_sigma <- quote(sigma[Mg] ~ dhalfnorm(0,1))
hiermod_out_dir <- "out/hiermod/16S_1_lognormal_MD"

fit_MDv <- readRDS(file.path(hiermod_out_dir, "fit_MDv.rds"))

## MDv -- Posterior contrast ---------------------------------------------------

pf_MDv <- post_full(fit_MDv, means_MDv)
pc_MDv <- compute_contrasts(pf_MDv, keep = "mean", group_levels = idx$Mg$levels)
save_report("fit_summary", model_id_MDv, fit_MDv, pc_MDv, model_MDv, model_name = "The Loose Cannon")

p_MDv_contrast <- contrast_plot_panels(pc_MDv, quant = c(0, 0.999), group_pal = Management_palette) +
  labs(x = 'Mean Hill number of order 1'); p_MDv_contrast

save_gg("fit_contrast_density", model_id_MDv, p_MDv_contrast, width = 8, height = 4)
