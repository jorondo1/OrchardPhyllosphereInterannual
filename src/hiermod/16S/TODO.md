# TODO list

Running list for the 16S lognormal hierarchical modeling work.

## Model 5 (MDLS2)

- [ ] Manuscript methods: short paragraph on the prior sensitivity check
  (5.5_MDLS2_16S_sensitivity.R) -- tight (dexp(2.5)) vs loose (dexp(1.5))
  sigma priors barely moved sigma[cell]/sigma_loc/sigma_tr or the May/July
  contrasts, so the calibrated prior isn't driving the result.
- [ ] Manuscript methods: short paragraph justifying sigma[cell] (residual
  SD split by Season, not just Management) -- see model comparison against
  the sigma[Mg]-only reduced fit in 5.6_MDLS2_16S_comparison.R.
