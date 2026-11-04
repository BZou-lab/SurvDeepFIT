# RSF feature-selection audit and diagnostic (2026-11-04)

Scope: Scenario 1, training n=1000, rho=0, ten repetitions. No complete
1,000-repetition study is authorized by this diagnostic. Cox, AFT, SurvDeepFIT,
and the existing simulation generator/configurations are unchanged.

## Implementation and provenance

The installed project library contains randomForestSRC 3.9.0 (R 4.5.0).
The audit inspected the installed function bodies, help, exports, and NEWS,
in addition to the reference documentation below. Installed source is decisive
where website documentation describes older interfaces or different defaults.

The v2 execution path is `run_simulation.R` -> `run_method()` in
`R/workflow.R` -> `fit_feature_selection()`/`fit_model()` in `R/fit_models.R`
-> `evaluate_features()` in `R/evaluate_features.R`. `summarize_results.R`
aggregates per-repetition rates and per-feature selection frequencies.

The original v2 RSF grow call is `rfsrc(Surv(time,status) ~ ., data, ntree=500)`.
It initially requests no importance. `subsample(fitted, B=100)` subsequently
computes importance with `importance=TRUE`, which resolves to `"anti"`.
Because importance is absent, subsample's default `block.size=1` is used for
the reference VIMP and inherited by its replicate forests. Subsampling refits
100 forests; these are not 100 permutations of a fixed fitted forest.

Effective survival defaults include mtry=ceil(sqrt(50))=8, nodesize=15,
logrank splitting, per-tree sampling without replacement of approximately
63.2% of available observations, and Uno-weighted concordance performance.
The subsample fraction is exp(-1) for n<=1000 and
exp(-1)*sqrt(1000/n) thereafter. At n=1000 each subsample has nominal size
367.9, with integer allocation within event/censoring strata. This fraction
is separate from each tree's sampling fraction. The informal description
"approximately inverse square root of n" omits the important scaling constant.

Configuration selects `extract.subsample(...)$var.sel$pvalue < .05`.
It does not select by positive importance alone. It does not use `signif`.

Scenario 1 has 50 candidate features: X1:X40 and Z1:Z10. The important set
is X1, X11, X21, X31, Z1; 45 others are null. The current generator uses
log-risk coefficients (0.1,0.2,0.3,0.4,0.5) with no sample-size scaling.
Those settings, Gompertz survival, Weibull censoring, and the independent
200-subject test sample are retained. Test outcomes do not enter selection.

Older scripts under `../numerical_experiments/` call `mod_permfit()` with
`method="random_forest", ntrees=500`, followed by `subsample()` and
`extract.subsample()`. `mod_permfit()` in
`../../../R/source_permfit_survival_20230610.r` forwards arguments to
`rfsrc()`. The old `ntrees` argument was ineffective; v2 already uses `ntree`.
Historical numerical-table aggregation/provenance is insufficient to identify
the exact p-value table and package version that produced the reported .74-.90.

The three saved v2 pilot RSF records found locally all report 3.9.0 and the
same single dataset/result: power .8, null-selection rate 16/45=.35556,
20 selected predictors, p-value rank sum 70. They are repeated executions,
not three independent simulation repetitions. Their elapsed times are
262.665, 282.300, and 382.145 seconds, each on one thread. Thus the historical
package bug below cannot by itself explain the inflation in these v2 pilots.

## Exact inference calculations in installed 3.9.0

For a feature let theta be full-data importance, theta_b its subsample
importance, n the training size, r the subsample fraction, and m=nr.

* `var.sel`: roots z_b=sqrt(m)*(theta_b-theta); empirical one-sided p-value
  mean(z_b > sqrt(n)*theta). There is no plus-one Monte Carlo adjustment;
  ties are excluded. With B=100 complete estimates, p-values lie on a .01
  grid and zero can occur. The strict p<.05 rule accepts 0,.01,.02,.03,.04.
  Nonparametric interval endpoints are theta minus reversed root quantiles
  divided by sqrt(n). Its `median` column is an interval median, not the
  original full-data VIMP.
* `var.sel.Z`: SE_sub=sqrt(r)*sd(theta_b); p=1-Phi(theta/SE_sub).
* `var.jk.sel.Z`: SE_jk=sqrt(r*mean((theta_b-theta)^2));
  p=1-Phi(theta/SE_jk). This includes the squared difference between the
  subsample mean and full-data estimate as well as subsample dispersion.

Both normal p-values are one-sided tests of positive VIMP. At extraction
alpha=.05 the displayed intervals are two-sided 95% intervals. Therefore
their lower>0 rule corresponds to one-sided p<.025, not p<.05. For a
one-sided .05 normal test the equivalent bound is theta-qnorm(.95)*SE>0;
the lower endpoint of the two-sided 90% normal interval is the same bound.
For nonparametric results, quantile interpolation and p-value discreteness
also cause boundary discrepancies. The diagnostic retains both decisions
for inspection without changing the prespecified p-value selection rule.

NEWS for 3.8.0 records a correction to nonparametric p-values to use the
centered subsampling distribution underlying the intervals, together with
alignment fixes. Installed 3.9.0 has that correction. Old results require
their original package version and saved table before claiming this was
the cause of inflated error.

The package uses `na.rm=TRUE` inside replicate summaries. The diagnostic
records available replicates and treats incomplete resampling or zero/nonfinite
standard errors as unavailable inference. It does not turn missing p-values
into nonselections or silently reduce null denominators. Finite zero empirical
p-values are retained when inference is otherwise estimable. Negative VIMP
is retained; it is not truncated or treated as an absolute importance value.

The existing feature evaluator correctly matches by predictor name and
computes selected signals/5 and selected nulls/45. All valid p-values in
each truth group are required for that group's rate. Partial `n_selected`
counts in the general evaluator can understate unknown selections; the
isolated diagnostic instead reports that count as unavailable if incomplete.
`setequal()` is not a duplicate check; the diagnostic adds an explicit
duplicate-row check. No mismatch or invalid p-value occurred in the saved pilots.

## Interpretation of alternative selection procedures

Permutation VIMP measures deterioration in predictive performance after
permuting a feature, whereas anti VIMP perturbs tree routing. Permutation
VIMP is a useful standard choice but its sign is not a significance test.
It does not inherently fix bias or establish type I error control.

An additional inferential distinction is relevant at rho=0. Under independence
of a null feature from both survival and the other predictors, independently
permuting that feature preserves their population joint distribution; for a
fixed predictor, its expected loss is therefore unchanged. Anti-split routing
does not preserve that distribution of routes. For example, flipping an
unbalanced null split can move many cases into a smaller, less accurately
estimated subtree. Thus positive anti importance need not identify a true
prognostic factor in a finite fitted forest. This is an inference from the
perturbation definitions, not proof that it explains the entire observed
inflation. It motivates permutation VIMP for this identification objective,
while leaving finite-sample OOB p-value calibration to be evaluated.

Subsampling inference is already present. Its validity relies on independent
subjects, a suitable limiting distribution/rate for the VIMP estimator,
and subsamples that approximate the sampling law. Normal inference adds a
normal approximation and adequate variance estimation. Survival inference
also requires appropriate censoring assumptions for the performance measure;
the simulation's independent censoring supports this. Finite forests,
training-size-dependent importance/bias, null degeneracy, and B=100 tail
instability can undermine calibration. The method is approximate, not an
exact conditional-independence or structural-feature-null test. The cited
paper motivates the jackknife bias adjustment, especially at low fractions;
it does not guarantee calibration for this particular experiment.

At correlated settings, a structural null covariate can be a predictive proxy.
Marginal VIMP and coefficient/conditional-null tests can target different
hypotheses. Even rho=0 does not automatically remove finite-forest bias.
Per-feature error is not family-wise error or false discovery rate. No
multiplicity correction or threshold calibration to these results is added.

`var.select()` was removed in 3.4.0 and is absent from installed 3.9.0.
Older `method="md"` selected by average first-split depth relative to a
model-derived noise-depth threshold. It was an algorithmic selection rule,
not a per-feature alpha=.05 test. `max.subtree()` remains available for
depth analysis, but introducing a new depth selector or a different package
would be a separately specified method. No package is installed/downgraded.

## Prespecified diagnostic and implementation

`configs/rsf_diagnostic.yaml` references the unchanged Scenario 1 configuration.
`run_rsf_diagnostic.R` reuses its original design grid, base seed 20261102,
L'Ecuyer-CMRG task streams, data substream 0, and RSF model/inference substreams
40/41. Each dataset is shared across all arms. Both fit types start from the
same model/inference streams; permutation draws can subsequently advance RNG
differently, so this is pairing by dataset, not a claim of identical internal
forests or subsample indices across importance types.

Primary arms, fixed before diagnostic results:

1. Original grow/subsample path, anti VIMP and `var.sel$pvalue < .05`.
2. Explicit permutation VIMP with block.size=1, same nonparametric rule.
3. The same permutation subsamples, `var.jk.sel.Z$pvalue < .05`.

The normal subsampling table `var.sel.Z` is a secondary sensitivity analysis
from the same saved estimates. No extra forest fitting is needed for either
normal rule. Standalone alternative runtimes include their shared permutation
fit cost; they must not be summed as if each required separate fitting.

The optional global-null diagnostic reuses covariates, censoring, split, and
event uniforms, replacing log-risk by zero using the inverse Gompertz
transformation. It has 50 nulls and no important features, so power and
signal-ranking metrics are NA. The main generator is not modified.

Outputs include per-repetition and feature-level results, full-data VIMP and
its ranking, pairwise signal/null ranking AUC (1 best, .5 chance), p-value
rank sums (smaller better), selected/null counts, valid replicate counts,
standard errors, interval endpoints, elapsed times, data, and provenance.
Monte Carlo SEs are calculated over independent dataset repetitions, rather
than treating correlated feature decisions as independent observations.
Ten repetitions diagnose gross problems; they cannot precisely validate .05.

Checkpoints preserve raw VIMP replicates and full-data VIMP, not every fitted
forest, limiting storage. Input/code/version/thread hashes prevent resuming
into incompatible outputs. Fit-stage warnings and failures are saved.
In the original post-fit VIMP path, the saved forest's `importance` metadata
can remain `"none"` even after VIMP has been calculated. The diagnostic's
`settings$importance` records that stored flag; the effective anti method
is determined by the inspected subsample/get.importance code and arm label.
An original copy of `R/fit_models.R` and the installed RSF function bodies
are also archived in the diagnostic output directory.

Only the RSF branch of `R/fit_models.R` changes: omitted `models$rsf$importance`
retains the original call; explicit `anti`, `permute`, or `random` requests
VIMP at grow time with block.size=1. Existing `selection$rsf_table` already
supports all three tables. No default production selection rule changes.

To request the provisional alternative in a separately copied production
configuration, add `models$rsf$importance: permute` and set
`selection$rsf_table: var.jk.sel.Z`. Keep alpha=.05. This is an explicit
scientific choice pending the diagnostic, not a blanket recommendation to
launch the full study.

Merge the following RSF settings into a copy of the complete study YAML;
retain the other keys in its `selection` and `models` mappings:

```yaml
selection:
  alpha: 0.05
  rsf_table: var.jk.sel.Z
models:
  rsf:
    ntree: 500
    subsamples: 100
    importance: permute
```

From this directory, deterministic validation (no forest fitting):

```bash
Rscript --vanilla validate_rsf_diagnostic.R
```

Sequential diagnostic, or one repetition for timing:

```bash
Rscript --vanilla run_rsf_diagnostic.R --config configs/rsf_diagnostic.yaml
Rscript --vanilla run_rsf_diagnostic.R --config configs/rsf_diagnostic.yaml --rep 1 --design signal
```

`slurm/run_rsf_diagnostic.slurm` selects the project R 4.5.0 library and
one thread per requested CPU. A Slurm array supplies the repetition through
`SLURM_ARRAY_TASK_ID`. Set an appropriate partition and time limit at submission;
do not run compute-heavy fits on the login node. `--design signal`,
`--design global_null`, or `--design both` restricts the diagnostic design.
`--output DIRECTORY` permits an independent run without overwriting existing
checkpoints. Completed fits are reused when resuming the same signature.

CSV summaries contain all completed records under the output directory.
Concurrent array tasks may refresh those CSVs in different orders; after all
tasks finish, refresh the final tables from the atomic RDS checkpoints without
refitting. From this directory, for the dated diagnostic executed in this audit:

```bash
Rscript --vanilla -e 'source("run_rsf_diagnostic.R"); out <- "results/rsf-diagnostic-20261104"; summarize_diagnostic(out, readRDS(file.path(out, "provenance.rds"))$signature)'
```

The dated run used the archived `executed_run_rsf_diagnostic.R`. Two array
tasks initially failed before fitting while concurrently replacing/reading
the shared provenance file; both were restarted with identical seeds/code.
The maintained runner now locks initialization and never rewrites an existing
manifest. A deterministic concurrent-initializer test covers this correction.
Source-path detection was also made compatible with `source()`, and the default
output is now the undated `results/rsf-diagnostic`. The four numerical helper
function bodies are identical to the archived execution version, verified by
comparison of their parsed/deparsed bodies. The maintained runner deliberately
has a different code hash; use a new output directory for a fresh rerun rather
than altering signatures on previously saved fits. Existing study code outside
the RSF branch is unaffected.

## References

* Ishwaran H, Lu M. Standard errors and confidence intervals for variable
  importance in random forest regression, classification, and survival.
  Statistics in Medicine. 2019;38:558-582. DOI:10.1002/sim.7803.
  https://www.ishwaran.org/papers/IL.StatMed.2019.pdf
* Package release history (also inspected in installed NEWS.md):
  https://cran.r-universe.dev/randomForestSRC/NEWS
* Subsampling reference (installed 3.9.0 help/function bodies determine
  precise defaults and formulas):
  https://www.randomforestsrc.org/reference/subsample.rfsrc.html
* Grow/VIMP definitions and survival performance:
  https://www.randomforestsrc.org/reference/rfsrc.html
* Minimal depth:
  https://www.randomforestsrc.org/articles/minidep.html
* Historical var.select interface (not available in installed 3.9.0):
  https://www.randomforestsrc.org/reference/var.select.rfsrc.html

## Diagnostic results and final recommendation

Results will be added after the approved compute jobs complete. No claim of
calibration or final method choice is made on the basis of the saved pilot.

Timing/reproduction check: the first diagnostic anti arm took 530.530 seconds
on compute node c0403 with one thread. Its generated dataset exactly matches
the saved v2 pilot. All five extracted objects (`se.Z`, `se.jk.Z`, `var.sel`,
`var.sel.Z`, `var.jk.sel.Z`) match the saved pilot exactly under `all.equal()`.
There are 100 finite subsample estimates for each of 50 features, no warnings,
and no invalid inference. This validates preservation of the original path.
For this one dataset, full-data VIMP ranking AUC is .94667, while the original
selector still chooses 16 nulls. Good ranking does not imply calibrated testing.
