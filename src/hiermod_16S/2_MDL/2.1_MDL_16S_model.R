# Idem to ITS model #2 (see for more doc)
# except model calibration

## Model specification ---------------------------------------------------------

source('src/hiermod_ITS/2_MDL/2.1_MDL_model.R')

# Model 1 + b[Lo]*sigma_loc
# Update ITS model with new priors 

model <- alist(
  likelihood   = Dv ~ dlnorm(mu, sigma[Mg]),
  linear_model = mu <- loga[Mg] + b[Lo]*sigma_loc,
  prior_loga   = loga[Mg] ~ dnorm(2,2),
  prior_b      = b[Lo] ~ dnorm(0,1),
  prior_sigma  = sigma[Mg] ~ dexp(1),
  prior_sigma_loc = sigma_loc ~ dexp(1)
)

# Mean backtransformation wrapper -----------------

# IDEM with ITS 'src/hiermod_ITS/1_MD/1.1_MD_model.R'

## Data-generating function ---------------------------------------------------

# IDEM with ITS 'src/hiermod_ITS/1_MD/1.1_MD_model.R'

# Prior simulator -------------------

# IDEM with ITS 'src/hiermod_ITS/1_MD/1.1_MD_model.R'