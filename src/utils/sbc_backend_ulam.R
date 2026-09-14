# adapters letting the SBC package's compute_SBC() fit rethinking::ulam() 
# models directly, as an alternative to this repo's
# own hand-rolled run_sbc()/save_sbc_report() (src/utils/sbc_helpers.R).
# See the evaluation this adapter was built to validate for the full
# cost/benefit reasoning.

SBC_backend_ulam <- function(model, ...){
  structure(list(model = model, args = list(...)), class = "SBC_backend_ulam")
}

SBC_fit.SBC_backend_ulam <- function(backend, generated, cores){
  do.call(rethinking::ulam,
          c(list(flist = backend$model, data = generated, cores = cores), backend$args))
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
