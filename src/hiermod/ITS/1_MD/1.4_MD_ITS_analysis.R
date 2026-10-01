# MODEL 1b (MDv), ITS: posterior contrast from the saved fit

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MDv_ITS, means_MDv()
model_MDv <- model_MDv_ITS
hiermod_out_dir <- "out/hiermod/ITS_1_lognormal_MD"

fit_MDv <- readRDS(file.path(hiermod_out_dir, "fit_MDv.rds"))

## MDv -- Posterior contrast ---------------------------------------------------

pf_MDv <- post_full(fit_MDv, means_MDv)
pc_MDv <- compute_contrasts(pf_MDv, keep = "mean", group_levels = idx$Mg$levels)
save_report("fit_summary", model_id_MDv, fit_MDv, model = model_MDv)
save_posterior_kable("results_report", model_id_MDv, pc_MDv)

p_MDv_contrast <- contrast_plot_panels(pc_MDv, quant = c(0, 0.999), group_pal = Management_palette) +
  labs(x = 'Mean Hill number of order 1'); p_MDv_contrast

save_gg("fit_contrast_density", model_id_MDv, p_MDv_contrast, width = 8, height = 4)
