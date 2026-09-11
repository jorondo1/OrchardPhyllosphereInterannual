# Code review — naming, dead code, Bayesian practice

One-time cleanup pass over `src/hiermod/` (2026-09-02). Applied fixes are
listed for the record; everything else is a suggestion, not yet done.

## Applied this pass

- **Rename**: `postcounts_MD_sim` -> `post_counts_MD_sim`
  (`1_ITS_lognormal_MD.R`) -- same script named the same function's output
  with/without an underscore.
- **Dead file**: `ITS_setup.R` moved to `archive/` -- sourced a
  `utils/hiermod_utils.R` that no longer exists (split into
  `hiermod_core.R`/`postcontrast_helpers.R`/`saver_functions.R`/
  `predictive_checks.R`/`sbc_helpers.R` since); zero references to it
  anywhere else in the repo, would error on `source()`.
- **De-duplication**: `cell_summary()` (identical 1-liner, copy-pasted in
  scripts 4, 5, 5-shifted) moved into `postcontrast_helpers.R`.
  `plot_ppc_season_contrast_stats()` (identical ~15-line block, same 3
  scripts) moved into `predictive_checks.R`; `patchwork` moved from 3
  separate in-script `library()` calls into `0_ITS_lognormal_SETUP.R`'s
  `pacman::p_load()`.
- **`save_fit()`** added to `saver_functions.R` (same pattern as
  `save_gg`/`save_report`) and called right after every real (non-
  simulation) model fit, scripts 1 through 5-shifted. Motivation: this
  session needed to independently sanity-check Model 5-shifted's contrast
  posteriors and found no `.rds` had ever been saved for any real fit --
  only `precis()`/plot outputs, meaning any later audit needs a full refit
  just to get raw draws back.
- **Comment trims**: the numbered-script "Model specification" narrative
  blocks for Models 3 and 4, and `models/model_5_MDLS2.R`'s header, cut from
  multi-paragraph prose to a few lines pointing at `MODEL_HISTORY.md` for
  the full why (that file now carries the complete reasoning trail, so the
  in-script versions were largely duplicates). `contrast_plot_panels()`'s
  doc comment (`postcontrast_helpers.R`) tightened the same way, gotchas
  kept.
- **Rename**: the model's own `s` -> `s_conv` (Season main effect at
  Conventional) and `c[Yr]` -> `yr[Yr]` (Year fixed effect), across
  `models/model_4_MDLS.R`'s `alist`/`means_MDLS()`/`sim_div_MDLS()`/
  `simulate_from_priors()`, `models/model_5_MDLS2.R`'s equivalents, and
  every call site in `4_ITS_lognormal_MDLS.R`/`5_ITS_lognormal_MDLS2.R`/
  `5_ITS_lognormal_MDLS2_shifted.R` (`check_recovery()` calls,
  `sim_div_MDLS2()` calls, `extracted_prior$s`). `yr` also brings Year in
  line with `tr`/`Tr`'s existing case-mirroring convention. Verified via a
  fresh compile+fit smoke test: `extract.samples()` returns `s_conv`/`yr`
  with no leftover `s`/`c`, and `means_MDLS2()`/`post_full()`/
  `simulate_from_priors()` all run clean against the renamed model. Pure
  rename (same priors, same structure) — no fit needs to be rerun for this
  alone, but any *new* real fit from here on will report `s_conv[...]`/
  `yr[...]` in `precis()` instead of `s[...]`/`c[...]`.

## Naming — suggestions, not applied

Core indices (`Mg`/`Lo`/`Tr`/`Cv`/`Mo`/`Yr`/`Dv`), pipeline objects
(`dat`/`idx`/`pf`/`pc_full`), and most model parameters (`loga`, `sigma`,
`b[Lo]`, `sigma_loc`, `gap_shift`, `tr[Tr]`, `sigma_tr`) are used
consistently everywhere. Left alone below because each would touch a model
parameter name or span many files — a bigger, riskier edit than this pass:

- **`md`** (scripts 4/5/5-shifted) is `pf$median`, easily misread as
  "model" given `model_MD`/`MODEL 1 MD` elsewhere in the same files.
- The parameter-name-vs-index-name convention is still mixed even after the
  `s`/`c` rename: `tr`/`Tr` and now `yr`/`Yr` mirror their index by case,
  but `b`/`Lo`, `loga`/`Mg`, `g`/`Lo` (Model 3) use an unrelated letter. Not
  wrong, just not uniform.
- `*_labels` objects live in different places per model (`mdl_labels` in
  the numbered script for Model 2, `mdlv_labels` in the model file for
  Model 3, no dedicated object at all for Models 4/5) — minor, low
  priority.

**Not renamed, not removed**: `sigma_to_cv()` (`hiermod_core.R`) has no live
callers — only `cv_to_sigma()` (its inverse) is used in any current script,
`sigma_to_cv()` only appears in `archive/`. Left in place rather than
deleted: it's the natural inverse of a function that *is* used, one line,
and plausibly useful later for reporting a fitted `sigma[cell]` back as an
interpretable CV.

## Bayesian terminology & practice

**No misuse found.** Grepped every `.R` file for "significant," "p-value,"
"confidence interval"/"CI," "null hypothesis," "statistically," and similar
— none appear as a claim about the science. The one real p-value in the
codebase (`sbc_helpers.R`'s `summarize_sbc()`, a KS test on SBC ranks) is
the standard Talts-et-al. calibration diagnostic, not a significance test.
`PI89`/`HPDI`/"True value"/"Observed" are used consistently and correctly
in every plot label and comment, including the ones built this session.

**Practice gaps found, not fixed here:**

- `save_report()` only ever captured `precis()` output (mean/sd/rhat/
  ess_bulk) — no divergent-transitions/treedepth/E-BFMI counts are logged
  anywhere for a real fit (confirmed directly: `fit_summary_MDLS2_shifted.txt`
  has none). `run_sbc()` already calls
  `fit@cstanfit$diagnostic_summary(diagnostics = "divergences")` for SBC
  replicates (`sbc_helpers.R:75`) — `save_report()` could do the same for
  the real fit. Not applied here since it's a real behavior addition, not a
  mechanical fix.
- `sigma[cell] ~ dexp(3)` (`models/model_5_MDLS2.R`'s `pr_sigma`) still
  carries an inline "revisit via prior-predictive check like model 4"
  comment. Model 5's Prior predictive check section exists and has run, but
  nothing records whether it specifically reconsidered this rate — worth
  either a one-line confirmation note or dropping the stale comment.
- `c[Yr]` staying an unpooled fixed effect (no `sigma_yr`) means the report
  can only ever give 3 year point-estimates, never a year-to-year variance
  — already tracked in `TODO.md` ("Pool `c[Yr]` eventually"), cross-referenced
  here for visibility.
- In the real Model 5-shifted fit, `sigma[4]` (Org-July) is markedly
  smaller (0.226) than the other 3 cells (0.70/0.665/1.05) and has lower
  ESS (~1,816 vs >10,000+ elsewhere, rhat 1.002 vs 1.000). One plausible
  explanation ties to the open "Organic-July bimodality" TODO item: a
  bimodal cell fit by a single lognormal component could show up as an
  artificially tight `sigma` with the extra spread partly absorbed into
  tree-level offsets instead — a hypothesis, not a confirmed diagnosis.
