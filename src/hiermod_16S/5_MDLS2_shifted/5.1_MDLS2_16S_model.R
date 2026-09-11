# 5.1_MDLS2_model.R 
# Based on ITS model #5
# Shifted, with management-season interaction ("cell")

source('src/hiermod_ITS/5_MDLS2/5.1_MDLS2_model.R')

## Model definition -------------------------------------------------

# IDEM to ITS Model 5 but we redefine it entirely,
# because it might need tuning
model <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[cell]), # Link function
  
  # Model:
  main_model = mu <- loga[Mg] + gamma*(Mo-1) + b[Lo]*sigma_loc + yr[Yr] + tr[Tr]*sigma_tr,
  gamma_def  = gamma <- s_conv + gap_shift*(Mg-1), # interactive gap

  # Old Priors:
  prior_loga = loga[Mg]   ~ dnorm(4,2), # May baseline per management, so we need to comfortably cover the data
  prior_b    = b[Lo]      ~ dnorm(0,1), 
  
  # New priors:
  prior_s    = s_conv     ~ dnorm(0,1), 
  prior_gs   = gap_shift  ~ dnorm(0,1), 
  prior_yr   = yr[Yr]     ~ dnorm(0,1), 
  prior_tr   = tr[Tr]     ~ dnorm(0,1), 
  
  # Hyperpriors:
  pr_sigma     = sigma[cell] ~ dexp(3), # Management-specific noise
  pr_sigma_loc = sigma_loc ~ dexp(2), 
  pr_sigma_tr  = sigma_tr  ~ dexp(2)  
)

# Backtransforming function --------------------------------------------------

# IDEM to ITS Model 5

## Data-generating function ---------------------------------------------------

# IDEM to ITS Model 5

# Prior simulator (SBC/prior-predictive glue) ---------------------------------

# IDEM to ITS Model 5
