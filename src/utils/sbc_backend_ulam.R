# adapters letting the SBC package's compute_SBC() fit rethinking::ulam() 
# models directly, as an alternative to this repo's
# own hand-rolled run_sbc()/save_sbc_report() (src/utils/sbc_helpers.R).
# See the evaluation this adapter was built to validate for the full
# cost/benefit reasoning.

SBC_backend_ulam <- function(model, ...){
  structure(list(model = model, args = list(...)), class = "SBC_backend_ulam")
}

# chains is always set to whatever `cores` compute_SBC() hands us for this
# fit, never taken from backend$args -- compute_SBC()'s own cores_per_fit
# defaults to 1 whenever n_sbc is large relative to available cores (true
# for any real SBC run, confirmed for both an 8-core laptop and a 48-core
# HPC node at n_sbc=100), so a fixed chains=2 in the backend would mean 2
# chains contending for 1 allocated core on every replicate -- oversub-
# scription, not a future::multisession vs mirai question. Pass
# cores_per_fit=N explicitly to compute_SBC() instead if you want every
# replicate to keep N>1 chains (at the cost of fewer replicates running
# concurrently) -- this backend will follow whatever it's given either way.
SBC_fit.SBC_backend_ulam <- function(backend, generated, cores){
  do.call(rethinking::ulam,
          c(list(flist = backend$model, data = generated, chains = cores, cores = cores), backend$args))
}

SBC_fit_to_draws_matrix.ulam <- function(fit){
  posterior::as_draws_matrix(fit@cstanfit$draws())
}

SBC_fit_to_diagnostics.ulam <- function(fit, fit_output, fit_messages, fit_warnings){
  # quiet = TRUE: this call itself would otherwise print the same
  # divergence/treedepth text a second time (on top of cmdstanr's own
  # automatic post-sample message) -- captured numerically below instead.
  ds <- fit@cstanfit$diagnostic_summary(diagnostics = c("divergences", "treedepth", "ebfmi"), quiet = TRUE)
  data.frame(n_divergent = sum(ds$num_divergent), n_max_treedepth = sum(ds$num_max_treedepth),
             n_low_ebfmi = sum(ds$ebfmi < 0.3))
}
