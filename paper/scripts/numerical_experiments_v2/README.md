# Numerical experiments, version 2

This framework implements the two manuscript scenarios for SurvDeepFIT, Cox,
Weibull AFT, RSF, and XGBoost. The original numerical scripts and the package's
`R/source_permfit_survival_20230610.r` are unchanged. No full simulations or SLURM
jobs have been run as part of this refactor.

## Files

```text
configs/scenario{1,2}.yaml    Study design and model/evaluation settings
R/generate_data.R            Shared Gaussian/Bernoulli features and survival outcomes
R/fit_models.R               Model adapters and method-specific feature selection
R/evaluate_features.R        Power, type I error, and ranking
R/evaluate_prediction.R      Held-out C-index, survival probabilities, and IPCW IBS
R/workflow.R                 Configuration, RNG streams, checkpoints, failure capture
run_simulation.R            One repetition or a local sequence of repetitions
summarize_results.R         Coverage, recovery, prediction, and failure tables
validate.R                  Small deterministic checks and inexpensive Cox/AFT fits
slurm/run_simulation.slurm   Job-array wrapper
results/                    Ignored generated outputs
```

## Environment and configuration

Use an existing R environment with `yaml`, `MASS`, and `survival`. RSF additionally
needs `randomForestSRC`; XGBoost needs `xgboost`. SurvDeepFIT uses `deepTL` and its
dependencies, plus `magrittr`, `Hmisc`, and `timeROC` for the unchanged legacy
inference function. Python, DeepSurv, DeepHit, SHAP, and JMIM are not needed for
these five methods. The scripts never install packages or call `check_dependency()`.
Use a compatible, validated `deepTL` build; its ensemble implementation is external
to this repository. Record/freeze that build before production.

The supplied YAML files differ only in `scenario`. They specify 1,000 repetitions,
training sizes 1,000/3,000, correlations 0/0.3/0.5, and block size 10 (50 predictors).
`test_fraction: 0.2` means **test size / training size**: 200/600 independent test
subjects, not a 20% split of total observations. Effects are scaled by
`1 / sqrt(n_train / 1000)`. The five/six important original features are recorded
in every dataset; generated interaction terms are not supplied as model inputs.

Copy a YAML file to configure a different study. `--n` and `--rho` select settings
already in its grid; `--repetitions` and `--seed` override those values. `--methods`
selects a comma-separated subset. `--rep` selects one repetition. Filtering methods
or settings does not change their assigned RNG streams. Both default scenario
configs share a grid and use disjoint stream IDs. If changing the design grid,
use a new output directory and treat it as a new study.

## Local checks and execution

From this directory:

```sh
# Configuration only; no data generation or model fitting.
Rscript run_simulation.R --config configs/scenario1.yaml --rep 1 --validate-only

# Small checks: generated fixtures, cheap Cox/AFT fits, metric and checkpoint tests.
# Writes ignored outputs under results/validation; does not train DNN/RSF/XGBoost.
Rscript validate.R

# Example after review: one manuscript-sized repetition, just Cox and AFT.
Rscript run_simulation.R --config configs/scenario1.yaml --rep 1 --n 1000 --rho 0 \
  --methods Cox,AFT --threads 1 --output results/study
```

For a smaller exploratory run, copy a config and reduce its training sizes, block
size, repetitions, and (if wanted) computational model settings. Such results are
smoke results, not manuscript estimates. Default SurvDeepFIT settings are expensive:
5 folds, 100 permutations, and ensembles of 100 DNNs. Running without `--rep`
processes every configured repetition; do not do that as a dependency check.

## Longleaf

Prepare the desired R environment using Longleaf's current module/environment
instructions. No module name, account, partition, GPU, or environment path is
hardcoded. Submit from this directory, or set `SIMULATION_ROOT` to its absolute
path. Use an absolute output path on appropriate persistent/scratch storage.

```sh
mkdir -p slurm/logs
sbatch --array=1-1000%20 --cpus-per-task=1 --mem=8G --time=12:00:00 \
  --output=slurm/logs/scenario1-%A_%a.out --error=slurm/logs/scenario1-%A_%a.err \
  slurm/run_simulation.slurm "$PWD/configs/scenario1.yaml" /path/to/study-results
```

The resource values above are examples, not measured requirements. Set memory,
walltime, CPUs, array concurrency, account, and partition through `sbatch` options.
Set `RSCRIPT` if an explicit R executable is needed. Add `--n 1000 --rho 0`, for
example, after the output argument to split the design into smaller submissions.
Submit scenario 2 separately with its config and distinct log filenames; it may
share the output root because scenario directories are separate. Set the array
range to the configured number of repetitions. Each array task handles one
repetition across its selected settings and methods. No GPU is requested.

Each task uses a L'Ecuyer-CMRG stream, with separate method/stage substreams.
Results record the actual RNG vectors, base seed, configuration, source hashes,
Git hash/dirty state, R/package versions, thread count, and session information.
Numerical reproducibility across different dependency builds or thread counts is
not promised. The legacy permutation function still resets its internal seeds;
the wrapper isolates those resets from other methods/stages but does not change
the inference routine's internal permutation coupling.

## Checkpoints and resuming

```text
results/study/scenario1/n1000/rho0/rep0001/
  data.rds                         One shared training/test dataset
  run.rds                          Configuration, RNG, and source provenance
  Cox/selection.rds                Fit, inference, features, and warnings
  Cox/prediction_full.rds           Prediction metrics and stage status
  Cox/prediction_selected.rds
  Cox/result.rds                   Compact, consistently structured result
  ... other methods ...
```

Rerun the same command to reuse successful stages and retry failed stages. A failed
prediction does not repeat successful feature inference. Selection/model objects
are checkpointed because regenerating ensembles or RSF inference can be expensive.
Convergence/numerical warnings reject a stage; other warnings are retained with
`status = warning`. Invalid feature p-values remain NA and are flagged; they are
never counted as successful non-rejections. C-index is retained if IBS alone fails.

A task-directory lock prevents simultaneous writes for the same setting/repetition.
On Longleaf checkpoints use same-directory temporary files and atomic rename.
Windows replacement keeps a `.previous` recovery file until replacement completes.
An interrupted task may leave `.lock`; verify that no task is running before
manually removing that empty lock directory and resubmitting. Do not delete data
or result files to force a rerun. Configuration, source, R/package versions, and
thread-count mismatches require a new output directory. Increasing only the number
of repetitions is allowed. Submit summaries after jobs finish, not during writes.

## Aggregation

```sh
Rscript summarize_results.R --input results/study --output results/study/summary
```

The script writes `summary.rds` and CSV data frames:

- `coverage`: configured repetitions, recorded results, missing repetitions, and
  selection failures for every scenario/size/correlation/method.
- `feature_recovery`: mean/SD power, per-null-feature type I error, and rank sum.
- `feature_rates`: feature-specific rejection/selection frequencies and ranks.
- `prediction_summary`: mean/SD C-index and IBS for full/selected models.
- `feature_results`, `prediction_results`, and `status`: repetition-level audit data,
  including warnings and failure messages.

Every aggregate reports its valid and missing counts; results with missing metrics
are not turned into zeros. Failed selections receive explicit missing rows.
XGBoost power/type I error remain NA because it supplies no p-values. A per-replicate
power or type I error is NA if any required feature p-value is invalid. Featurewise
rates still report their own denominators. Missing repetitions are separate from
failed recorded attempts. Do not present an incomplete subset as 1,000 repetitions.
Duplicate method/repetition records or mixed signatures within a scenario stop
aggregation. Tables are data frames/CSV, allowing final manuscript formatting
without manual recalculation.

## Equivalence, corrections, and decisions before production

The distribution, feature effects, effect scaling, architecture, alpha=0.05, and
requested model settings follow the manuscript. This is distributional equivalence,
not reproduction of historical random numbers or previously published results.

Confirmed implementation problems addressed locally in v2:

- Old censoring generation rewound the covariate RNG. The shared generator now
  draws covariates, event times, and censoring sequentially. Exact inverse-CDF
  Gompertz sampling replaces `simsurv`, avoiding its internal seed reset.
- Old prediction scripts read PermFIT outputs for every method, despite manuscript
  descriptions of native selection. Cox/AFT now use native coefficient p-values;
  RSF uses subsampling; XGBoost selects the ten largest gain scores. AFT predictor
  rows are matched by name, excluding intercept and log(scale). These extraction
  rules do not establish what the missing historical aggregation code did.
- RSF now receives `ntree`, not the ineffective `ntrees`. Its current 500-tree
  default may have masked that typo. XGBoost explicitly receives 50 rounds,
  depth 3, eta 0.10, and the Cox objective from Supplementary Note 3.
- AFT discrimination uses negative predicted log-time (the original direction
  was already correct). Weibull probabilities use the fitted AFT scale. Weibull
  misspecification under Gompertz survival is retained as a methodological choice;
  no robust covariance, alternative AFT family, or new inference is introduced.
- New local Breslow prediction uses the inclusive risk set `time >= event_time`,
  correcting the package helper's strict `>` denominator. The package helper itself
  is unchanged. Native RSF survival curves and native Weibull AFT probabilities
  are used instead of treating arbitrary risk scores as Cox linear predictors.
- Failed fits/warnings are checkpointed explicitly, not saved as apparently
  successful `try-error` objects. No invalid-prediction observations are dropped
  silently from metrics. Test outcomes are used only for evaluation/IPCW weights.

The following need scientific review before a full simulation launch:

1. **SurvDeepFIT inference remains legacy and unvalidated.** It calls the unchanged
   core `permfit_survival()` and selects `importance_pval_x`, matching the original
   prediction scripts. Its permutation-variability denominator differs from
   manuscript Algorithm 1 and does not establish a population sampling standard
   error. No variance estimator, p-value, fold procedure, or internal permutation
   seed was changed. Reconcile the manuscript and inference code before interpreting
   power/type I error. External `deepTL` bootstrap/filtering also needs version review.
2. **RSF inference table:** `var.sel` matches the companion native-method application
   code. Supplementary Note 2 does not identify the exact subsampling table; confirm
   the historical table and package version. `subsamples=100` makes the current
   default explicit. Other RSF tuning settings remain package defaults.
3. **No sparse-selection fallback:** numerical manuscript text says p<0.05. Old
   scripts instead replaced fewer than two selected variables with the top six.
   V2 does not add nonsignificant variables. With zero selected variables the selected
   prediction is explicitly unavailable, while the full model is still evaluated.
   Confirm the intended policy before comparison with old prediction figures.
4. **IBS details were not specified and numerical IBS producers were missing.** V2
   exposes a common 100-point grid between the training event-time 10th/90th
   percentiles, test-set censoring KM, IPCW with G(T-) for events, and trapezoidal
   integration divided by interval length. Unsupported censoring weights give NA,
   not arbitrary clipping. Review the interval against the historical IBS analysis.
5. **Ranking ties:** p-value methods rank ascending p-value; XGBoost ranks descending
   gain, with unused features assigned zero gain. Ranks use average ties, and top-ten
   gain ties break by feature name. The manuscript does not specify tie handling.
   Per-null-feature type I error is not family-wise error; no multiplicity correction
   was added. Marginal permutation importance is not automatically a conditional
   structural-null test when predictors are correlated.

The refactor cannot establish why historical AFT/RSF type I errors were extreme:
the numerical table aggregation code and one Scenario 2 generator were absent.
Both scenarios/settings are now supported by the common generator, but original
results still need provenance checks. Do not change competitors merely to improve
their results.

## Validation scope

`validate.R` checks both mechanisms at each correlation on small fixtures, exact
outcome construction, reproducible/disjoint stage streams, feature metrics with
missing p-values, C-index direction, an exact inclusive-risk-set Breslow fixture,
uncensored IBS, inexpensive native Cox/AFT inference/predictions, warning rejection,
failed-stage retry, successful-stage reuse, signature protection, and save/load/CSV
aggregation. It does not statistically validate the legacy p-values, estimate power,
or run ensembles/forests/boosting. Validate those package-specific adapters in the
intended Longleaf environment before running the full study.

Local verification on 2026-11-02 passed under R 4.4.0: the validation script,
configuration-only CLI, and an end-to-end 120-training/24-test-subject Cox/AFT run
followed by CLI aggregation. The working local environment could not load
SurvDeepFIT's dependencies, RSF, or XGBoost; those adapters were not trained.
The Windows sandbox also prevented Bash from starting for a syntax-only check,
so the SLURM wrapper still needs checking with Bash on Longleaf. No packages were
installed and no jobs were submitted.
