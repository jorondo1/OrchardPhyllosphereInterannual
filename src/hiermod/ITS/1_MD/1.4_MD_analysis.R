# MODEL 1 (MD) + MODEL 2 (MDv): posterior contrast, run against the saved fits.

hiermod_marker <- "ITS"
source('src/hiermod/0_SETUP.R')
source('src/hiermod/Models/MD_model.R') # model_MD_ITS/model_MDv_ITS, means_MD/means_MDv
model_MD  <- model_MD_ITS
model_MDv <- model_MDv_ITS
hiermod_out_dir <- "out/hiermod/ITS_1_lognormal_MD"

fit_MD  <- readRDS(file.path(hiermod_out_dir, "fit_MD.rds"))
fit_MDv <- readRDS(file.path(hiermod_out_dir, "fit_MDv.rds"))

## MD -- Posterior contrast ---------------------------------------------------
# Back-transform loga/sigma posterior draws to a raw-scale group contrast

pf_MD <- post_full(fit_MD, means_MD)
pc_MD <- compute_contrasts(pf_MD, keep = "mean", group_levels = idx$Mg$levels)
save_report("fit_summary", "MD", fit_MD, pc_MD, model_MD, model_name = "The Bare Bones")

p_MD_contrast <- contrast_plot_panels(pc_MD, quant = c(0, 1), group_pal = Management_palette) +
  labs(x = 'Mean Hill number of order 1'); p_MD_contrast

save_gg("fit_contrast_density", "MD", p_MD_contrast, width = 8, height = 4)

## MDv -- Posterior contrast ---------------------------------------------------

pf_MDv <- post_full(fit_MDv, means_MDv)
pc_MDv <- compute_contrasts(pf_MDv, keep = "mean", group_levels = idx$Mg$levels)
save_report("fit_summary", "MDv", fit_MDv, pc_MDv, model_MDv, model_name = "The Loose Cannon")

p_MDv_contrast <- contrast_plot_panels(pc_MDv, quant = c(0, 0.999), group_pal = Management_palette) +
  labs(x = 'Mean Hill number of order 1'); p_MDv_contrast

save_gg("fit_contrast_density", "MDv", p_MDv_contrast, width = 8, height = 4)
