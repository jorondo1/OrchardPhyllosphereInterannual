# 1.3_MD_analysis.R -- MODEL 1 (MD) + MODEL 2 (MDv): posterior contrast,
# run against the fits saved by 1.2_MD_validation.R -- no refit needed.

source('~/Repos/orchardPhyllosphere2/src/hiermod/ITS/0_SETUP.R')
source('~/Repos/orchardPhyllosphere2/src/hiermod/ITS/1.1_MD_model.R') # model_MD/model_MDv, means_MD/means_MDv, postcounts_Model1/2()
hiermod_out_dir <- "out/hiermod/ITS_1_lognormal_MD"

fit_MD  <- readRDS(file.path(hiermod_out_dir, "fit_MD.rds"))
fit_MDv <- readRDS(file.path(hiermod_out_dir, "fit_MDv.rds"))

dat_MD <- list(
  Dv = div$Hill_1,
  Mg = idx$Mg$to_index(div$Management)
)

## MD -- Posterior contrast ---------------------------------------------------
# Back-transform loga/sigma posterior draws to a raw-scale group contrast

post_counts_MD <- postcounts_Model1(fit_MD)
save_report("fit_summary", "MD", fit_MD, post_counts_MD, model_MD, model_name = "The Bare Bones")

colnames(post_counts_MD) <- c(idx$Mg$levels, 'Contrast')

p_MD_contrast <- plot_contrast_density(post_counts_MD, group_name = 'Posterior mean') +
  labs(x = 'Mean Hill number of order 1') +
  scale_colour_manual(values = Management_palette) +
  scale_fill_manual(values = Management_palette) +
  theme(legend.position = c(0.8, 0.8)); p_MD_contrast

save_gg("fit_contrast_density", "MD", p_MD_contrast, width = 8, height = 4)

## MDv -- Posterior contrast ---------------------------------------------------

post_counts_MDv <- postcounts_Model2(fit_MDv)
save_report("fit_summary", "MDv", fit_MDv, post_counts_MDv, model_MDv, model_name = "The Loose Cannon")

colnames(post_counts_MDv) <- c(idx$Mg$levels, 'Contrast')
p_MDv_contrast <- plot_contrast_density(post_counts_MDv, quant = 0.999, group_name = 'Posterior mean')+
  labs(x = 'Mean Hill number of order 1') +
  theme(legend.position = c(0.9, 0.8)); p_MDv_contrast

save_gg("fit_contrast_density", "MDv", p_MDv_contrast, width = 8, height = 4)
