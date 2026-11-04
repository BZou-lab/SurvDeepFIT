# Numerical experiments, version 2

The approved, isolated RSF audit/diagnostic is documented in
[`RSF_AUDIT.md`](RSF_AUDIT.md). Run `run_rsf_diagnostic.R` with
`configs/rsf_diagnostic.yaml` to compare the original anti/nonparametric
selector with permutation/nonparametric and permutation/jackknife selectors
on paired datasets. Both study configurations now use permutation importance
and nonparametric subsampling p-values (`var.sel`, p < 0.05). The diagnostic
explicitly preserves the original anti arm for comparison.

This framework implements the two scenario structures for SurvDeepFIT, Cox,
Weibull AFT, RSF, and XGBoost. The current generator uses the manuscript effect
coefficients with no sample-size scaling, doubled from the previous pilot design.
The original numerical scripts and the package's
`R/source_permfit_survival_20230610.r` are unchanged. A single-repetition Longleaf
pilot has been reported; the full simulation study has not been launched here.

## Files

```text
configs/scenario{1,2}.yaml    Study design and model/evaluation settings
R/generate_data.R            Shared Gaussian/Bernoulli features and survival outcomes
R/fit_models.R               Model adapters and method-specific feature selection
R/evaluate_features.R        Power, type I error, and ranking
R/evaluate_prediction.R      Held-out C-index, survival probabilities, and IPCW IBS
R/workflow.R                 Configuration, RNG streams, checkpoints, failure capture
run_simulation.R            One repetition or a local sequence of repetitions
check_environment.R         Dependency loading and interface checks without fitting
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
subjects, not a 20% split of total observations. The current effects are:

- Scenario 1: `0.2*X1 + 0.4*X11 + 0.6*X21 + 0.8*X31 + Z1`.
- Scenario 2: `Z1 + 2*X1*Z2 + 0.5*X11^2 + X21*X31`.

The former `1 / sqrt(n_train / 1000)` scaling is disabled, so both training
sizes use the same coefficients. The absence of scaling differs from the manuscript.
Existing pilot and dated RSF diagnostic outputs used the previous, halved effects.
The five/six important original features are recorded
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

### Single repetition pilot

Edit locally, review your diff, and commit/push the intended changes to GitHub on
`simulation-refactor`. Uncommitted local edits cannot be pulled on Longleaf. In
the Longleaf clone, first check for local changes, then update the same branch:

```sh
git status
git fetch origin
git switch simulation-refactor
git pull --ff-only origin simulation-refactor
git rev-parse HEAD
cd paper/scripts/numerical_experiments_v2
```

If `git status` shows local edits on Longleaf, resolve them before pulling; do not
discard them automatically. Git transfers code, not your R environment. Load the
R module/environment that matches your installed packages and confirm `Rscript
--version`. Find available modules using `module spider r` if needed. The wrapper
uses that environment, or the explicit executable in `RSCRIPT`; it does not select
or install an R version. Repository attributes keep R/YAML/SLURM text files in LF
format when Git checks them out on Linux.

Check configuration and dependencies without simulation or model fitting:

```sh
Rscript run_simulation.R --config configs/scenario1.yaml --rep 1 --n 1000 --rho 0 --validate-only
Rscript check_environment.R --config configs/scenario1.yaml
```

The dependency check prints package versions, load failures, and library paths.
For SurvDeepFIT, it also runs hand-built `dnnet`/`dnnetEnsemble` prediction fixtures
through the real deepTL methods, legacy prediction helper, and v2 prediction adapter.
It additionally exercises the real legacy AUC/IPCW survival-formula evaluation path
on a fixed censored fixture, checking AUC and all three inference score differences.
It performs no model training or simulated-data generation. To check only native Cox/AFT requirements, add
`--methods Cox,AFT`. The SLURM wrapper repeats the dependency check on its compute
node and stops before fitting if it fails.

Submit one repetition at one setting, using a new results directory for this
previous halved-effect design. Replace the example output path below before running:

```sh
export SIM_RESULTS="/path/to/persistent/survdeepfit-v2-pilot"
mkdir -p slurm/logs
sbatch --array=1 --cpus-per-task=1 --mem=8G --time=12:00:00 \
  --output=slurm/logs/pilot-%A_%a.out --error=slurm/logs/pilot-%A_%a.err \
  slurm/run_simulation.slurm "$PWD/configs/scenario1.yaml" "$SIM_RESULTS" \
  --n 1000 --rho 0
```

This runs repetition 1 of Scenario 1, with 1,000 training/200 test subjects,
correlation zero, and all five methods at their configured settings. Append
`--methods Cox,AFT` for a cheaper first execution check. Use `scenario2.yaml` for
the nonlinear pilot. Without `--n` and `--rho`, even a one-task array runs all
six settings. The pilot does not estimate power or type I error reliably.

Monitor the job ID returned by `sbatch` and inspect resource use after completion:

```sh
squeue -u "$USER"
sacct -j JOB_ID --format=JobID,State,Elapsed,MaxRSS,ExitCode
tail -n 80 slurm/logs/pilot-JOB_ID_1.out
tail -n 80 slurm/logs/pilot-JOB_ID_1.err
Rscript summarize_results.R --input "$SIM_RESULTS" --output "$SIM_RESULTS/summary"
```

Replace `JOB_ID` with the numeric job ID. Check `summary/status.csv` for selection
failures/warnings, `summary/prediction_results.csv` for prediction/IBS failures,
and `summary/coverage.csv` for missing methods/repetitions. A SLURM exit code of
zero does not guarantee every method succeeded: fitting failures are recorded so
other methods can finish. A pilot will intentionally show missing repetitions
against the configured 1,000-repetition study. Rerunning the identical submission
retries failed stages and reuses successful ones. Keep its source, settings,
environment, and thread count unchanged while it runs; use a new output directory
after code changes. Size later resource requests from the pilot's measurements.

### Correction after pilot job 8066047

The reported second pilot passed Cox, AFT, RSF, and XGBoost, but SurvDeepFIT selection
failed with `could not find function "Surv"`. The earlier S4 prediction fix is
preserved. A reduced local training run reproduced this error. An error calling
handler captured the live stack before `capture_stage()` could unwind it:

```text
permfit_survival -> log_lik_diff -> timedAUC -> timeROC
  -> pec::ipcw -> ipcw.marginal -> prodlim::prodlim
  -> EventHistory.frame -> model.frame.default -> eval
Failing call: Surv(failure_time, status)
```

In timeROC 0.4.1, `timeROC()` constructs `Surv(failure_time, status) ~ 1`
inside its own execution frame. Formula evaluation follows that frame's lexical
environment, the timeROC namespace, which cannot resolve `Surv` in the local build.
It does not follow the caller's legacy environment, so the existing `core$Surv`
binding cannot fix this path. The dependency's formula callback fails even though
network training, S4 prediction dispatch, and legacy C-index evaluation work.

The v2 adapter now binds a local copy of `timeROC()` whose environment is a private
child of the original function environment, containing `survival::Surv`. Its body
and formals are unchanged; existing package helpers remain accessible through the
original environment. The installed namespace and global search path are untouched.
No global package attachment is required. This repairs name lookup without changing
AUC/IPCW calculations, inference mathematics, model settings, data generation, or RNG
behavior. AFT and the other successful methods are unchanged.

Both `check_environment.R` and `validate.R` now execute the actual failing
`timedAUC`/IPCW/formula path on 24 fixed subjects with censoring, checking perfect
and reversed rankings (AUC 1 and 0), and all three `log_lik_diff` components.
Validation also checks that the package function environment/body/formals, global
search path, and RNG state remain unchanged by this repair and fixture.

Local verification on 2026-11-03 passed under R 4.6.1, survival 3.8.6,
deepTL 1.0.0, timeROC 0.4.1, Hmisc 5.3.0, and magrittr 2.0.5:
`validate.R`, the SurvDeepFIT-only preflight, and a real end-to-end SurvDeepFIT run
followed by CLI aggregation. The smoke run used 120 training/24 test subjects,
block size 2 (10 predictors), rho 0, 2 folds, 4 permutations, 3 networks per
ensemble, one hidden layer of 4 nodes, 10 epochs, batch size 10, and early-stop
setting 10. IBS used 20 points over training event quantiles 0.1 to 0.6.
Selection and both prediction/IBS stages were `ok`, with no stage warnings;
full/selected C-index was 0.8077/0.5000 and IBS was 0.12135/0.12174.
These are execution checks, not statistical validation of inference or power.
Local Windows startup emitted locale warnings before script execution.
Longleaf R 4.5.0/survival 3.8.3 verification remains necessary. No packages were
installed and no cluster jobs or full study were launched here.

After pushing the reviewed local changes and pulling them on Longleaf, run the
following from the outer workspace root. Keep the old pilot directories intact:

```sh
source scripts/seer_r_environment.sh
cd SurvDeepFIT_repo/paper/scripts/numerical_experiments_v2
Rscript --vanilla check_environment.R --config configs/scenario1.yaml --methods SurvDeepFIT
Rscript --vanilla validate.R
mkdir -p results slurm/logs
PILOT_OUTPUT="$(mktemp -d "$PWD/results/pilot-surv-scope-XXXXXXXX")"
echo "$PILOT_OUTPUT"
sbatch --array=1 --cpus-per-task=1 --mem=8G --time=12:00:00 \
  --output=slurm/logs/pilot-surv-scope-%A_%a.out \
  --error=slurm/logs/pilot-surv-scope-%A_%a.err \
  slurm/run_simulation.slurm "$PWD/configs/scenario1.yaml" \
  "$PILOT_OUTPUT" --n 1000 --rho 0 --methods SurvDeepFIT
```

The array index supplies repetition 1; the unchanged config supplies seed 20261102
and 200 test subjects. Only SurvDeepFIT runs, at the configured pilot settings.
Use a new directory because changed code invalidates old checkpoint signatures.
After completion, summarize with `Rscript --vanilla summarize_results.R --input
"$PILOT_OUTPUT" --output "$PILOT_OUTPUT/summary"` and inspect selection, prediction,
and IBS status; a zero SLURM exit code alone does not establish success.

### Corrections after pilot job 8039199

The reported pilot passed Cox, RSF, and XGBoost, but failed SurvDeepFIT prediction
dispatch and AFT convergence. The fixes preserve the generator, fitted families,
hyperparameters, feature-selection rules, and legacy SurvDeepFIT inference:

- **SurvDeepFIT:** namespace-only loading left unqualified legacy `predict()` calls
  resolving to the S3 generic. V2 now binds the S4 prediction generic retrieved from
  the deepTL namespace in the sourced legacy environment, and uses the same generic
  for its own risk/survival predictions. The preflight fixtures catch this failure
  without training a network or attaching packages to the global search path.
- **AFT:** the full Weibull model now starts from a checked intercept-only fit,
  with predictor coefficients initialized to zero and the fitted log-scale as its
  starting scale parameter. Both fits must pass the existing warning guard and
  finite-parameter/covariance checks. The likelihood, distribution, covariance
  estimator, default 30-iteration budget, and convergence tolerances are unchanged.
  Selection and all prediction refits use this helper. Increasing iterations alone
  or accepting a nonconverged model is not the fix.
- **Diagnostics:** AFT prediction rows now record initialization, iteration count,
  and maximum absolute likelihood score. Full initialization/control details are
  retained in the selection model's `aft_solver` attribute. RSF still uses native
  subsampling and `ntree`; the returned tree count is now checked against the request.
- **Environment:** the SLURM wrapper uses `Rscript --vanilla` for both preflight and
  simulation. Activate the intended R environment first, with its library paths
  exported (for example `R_LIBS_USER`). User/site R profiles and Renviron files are
  ignored so an old `.Rprofile` cannot replace those paths. No user startup files
  are edited. Use `Rscript --vanilla` for manual cluster checks too.

Local reconstruction of the reported design/seed reproduced the rounded Cox
C-index (0.6641). The default AFT initialization failed at 30 iterations; the new
initialization converged in 6, with maximum score below 10^-6. A 100-iteration
verification gave the same coefficients. Full/selected AFT prediction and IBS
completed successfully locally. This reconstruction is not a check of the remote
saved RDS, which is not present locally; Longleaf still needs the corrected pilot.
The earlier AFT/RSF review note's initialization concern is therefore supported.
Its separate legacy AFT `...` forwarding issue remains in the original package
helper, which v2 native AFT fitting does not call. No package-level patch was made.

After pushing locally and pulling on Longleaf, activate the same project R 4.5
environment used for job 8039199. From this v2 directory, repeat only the two
previously failed methods, in a **new** output root:

```sh
Rscript --vanilla check_environment.R --config configs/scenario1.yaml --methods SurvDeepFIT,AFT
Rscript --vanilla validate.R
mkdir -p results slurm/logs
PILOT_OUTPUT="$(mktemp -d "$PWD/results/pilot-fixed-XXXXXXXX")"
export PILOT_OUTPUT
sbatch --array=1 --cpus-per-task=1 --mem=8G --time=12:00:00 \
  --output=slurm/logs/pilot-fixed-%A_%a.out --error=slurm/logs/pilot-fixed-%A_%a.err \
  slurm/run_simulation.slurm "$PWD/configs/scenario1.yaml" "$PILOT_OUTPUT" \
  --n 1000 --rho 0 --methods SurvDeepFIT,AFT
```

Omit `--methods SurvDeepFIT,AFT` to repeat all five methods. The new root recreates
the shared dataset with the same design/seed. Keep the original pilot for reference;
do not aggregate old and corrected roots together. Unrequested methods will be
listed as missing in coverage. After completion, use the summary command above
with `$PILOT_OUTPUT`, and inspect selection, full/selected prediction, and IBS
statuses; do not infer success from SLURM's exit code alone. Do not launch the full
1,000-repetition study on the strength of this one pilot.

### Full study after pilot review

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

The feature distributions, scenario term structure, architecture, alpha=0.05,
and requested model settings follow the manuscript. The current manually edited
generator uses the original effect coefficients and removes sample-size scaling. Thus
its event-time distribution intentionally differs from the original manuscript
design. Validation checks these coefficients independently against
the original formulas, without the old scaling. Historical random numbers
and published results are not reproduced.

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
  The pilot correction changes numerical initialization only, as described above.
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
   default explicit. Both scenarios explicitly request `importance: permute`,
   using `block.size=1` and nonparametric `var.sel` p-values < 0.05.
   Other RSF tuning settings remain package defaults.
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

The native Cox/AFT adapter checks use a separate deterministic Weibull fixture,
and AFT is additionally checked on the reconstructed n=1,000, rho=0 pilot dataset.
The revised halved/unscaled Gompertz smoke dataset produced an AFT convergence
warning locally on 2026-11-03, and the failure guard correctly rejected it.
The model and convergence policy were not changed to accommodate that failure.
The scenario-generator checks still verify the actual revised Gompertz mechanisms.

Local verification on 2026-11-02 passed under R 4.4.0: the validation script,
configuration-only CLI, and an end-to-end 120-training/24-test-subject Cox/AFT run
followed by CLI aggregation. The working local environment could not load
SurvDeepFIT's dependencies, RSF, or XGBoost; those adapters were not trained.
On 2026-11-03, the updated validation suite and native-only dependency check passed,
and Bash's syntax-only check passed for the SLURM wrapper. The all-method dependency
check correctly rejected incompatible local compiled dependencies before fitting.
The revised suite passed with survival 3.5.8 under local R 4.4 and survival 3.8.6
under R 4.6. The R 4.6 environment also passed real deepTL prediction fixtures and
all-method preflight (deepTL 1.0.0, Hmisc 5.3.0, timeROC 0.4.1, randomForestSRC 3.7.0,
xgboost 3.2.1.1). No ensemble, forest, or boosting training was used for these fixes.
Corrected execution on Longleaf remains to be tested. No packages were installed
and no jobs were submitted here.
