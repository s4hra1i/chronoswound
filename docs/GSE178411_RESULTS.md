# GSE178411 prospective analysis results

## Provenance and question

This analysis tests whether a locked 15-gene transcript panel provides practically useful
resolution of time since injury among surgically sampled burn wounds collected 3–27 days after
injury. The protocol was published as `gse178411-protocol-v1.0` before any expression–outcome
model was examined. The analysis used implementation commit `b56ae32`; results were generated
with Python 3.12.13 and scikit-learn 1.8.0.

The claim is deliberately narrow. This is an early- versus late-subacute burn-wound analysis,
not general wound dating and not evidence for wounds sampled before 72 hours, post-mortem
injuries, other mechanisms or forensic casework.

## Cohort and evaluation

The identical four-model comparison contained 49 samples from 39 patients. Samples span 3–27
days after injury. Five-fold patient-grouped outer validation was repeated across 20 fixed seeds;
ridge tuning occurred only inside each outer training set. All 100 fold results, 20 pooled
repetition results, 3,920 sample-model-repetition predictions and patient-to-fold assignments are
retained in the results directory.

The primary comparison was a locked-panel ridge against a median calculated only from each
outer training fold. The prospective rule required all of: at least 20% lower mean repeated-CV
MAE; a patient-clustered BCa 95% interval for the paired improvement wholly above zero; and a
one-sided full-pipeline patient-block permutation *p* below 0.05.

## Results

| Model | Mean pooled MAE | Median | IQR | Range across 20 seeds |
|---|---:|---:|---:|---:|
| Training-fold median | 4.26 d | 4.23 d | 4.08–4.35 d | 4.08–4.63 d |
| Covariates only | 4.58 d | 4.55 d | 4.48–4.67 d | 4.30–5.05 d |
| Locked panel | **2.80 d** | 2.79 d | 2.65–2.90 d | 2.54–3.30 d |
| Covariates + panel | 3.03 d | 3.01 d | 2.92–3.09 d | 2.60–3.46 d |

The panel reduced MAE by 34.3% relative to the training-fold median. The paired improvement was
1.46 days, with a patient-clustered BCa 95% interval of 0.68–2.43 days. None of 199 complete
patient-block permutation reruns produced an improvement at least this large, giving the minimum
attainable one-sided Monte Carlo *p*=0.005. The three-part prospective success rule was therefore
met.

The preregistered covariates-only versus combined contrast favoured the combined model by 1.55
days (patient-clustered BCa 95% CI 0.88–2.24). This contrast is weak evidence of incremental
molecular value because the covariates-only model itself underperformed the training-fold median.
The more informative post-hoc comparison went in the opposite direction: adding the recorded
covariates worsened panel MAE from 2.80 to 3.03 days. The data therefore support signal in the
locked molecular panel, but provide no evidence that these clinical covariates improve its
temporal estimates.

The panel's descriptive Early/Late threshold accuracy was 75.5% and balanced accuracy 75.1%,
compared with 55.1% accuracy from always predicting the majority class.
These values have no success criterion and cannot strengthen or replace the regression result,
because the classes are a deterministic threshold of the continuous target.

## Interpretation

The result clears a deliberately conservative, prospectively fixed threshold and is stable
across the 20 registered fold seeds. It is substantially more informative than the earlier
GSE8056 analysis because GSE178411 contains individual samples, exact days and patient IDs.
Nevertheless, it remains internal resampling of one small cohort rather than external
validation.

Several constraints limit the inference. The 49 observations come from only 39 patients, and
31 patients contribute a single eligible sample. The target has 16 distinct values and is
concentrated around operation scheduling, so apparent temporal signal may partly reflect the
clinical process determining when tissue was collected. Sample-level collection year and RNA
integrity are unavailable despite collection spanning 2002–2018; storage-associated degradation
could therefore masquerade as temporal structure. Burn severity, treatment and other unrecorded
clinical factors may also confound the result.

Inner-CV ridge alphas were chosen by giving each inner fold equal weight rather than weighting
fold MAE by fold size. In the patient-clustered bootstrap, the estimand is the mean sample-level
error across repeated-CV rows: resampled patients retain all their observations, so patients with
multiple eligible samples contribute proportionally more rows. Both choices are retained as part
of the prospectively executed analysis and are stated here rather than changed after observing
the result.

The correct conclusion is that this locked panel demonstrated useful internal temporal
resolution for 3–27-day surgical burn specimens under the registered evaluation. It is not a
validated wound-age estimator and must not be used in real clinical, legal or forensic cases.

## Reproducible outputs

- `analysis_summary.json`: protocol, code, software, settings, model metrics and decisions;
- `oof_predictions.csv`: every sample-level out-of-fold prediction;
- `repetition_metrics.csv`: pooled MAE for each model and seed;
- `fold_metrics.csv`: all 100 outer folds for each model;
- `fold_assignments.csv`: patient-to-fold allocation for every repetition;
- `permutation_null.csv`: all 199 full-pipeline null improvements;
- `analysis_cohort.csv`: accessions, patient groups and eligible outcomes;
- `figures/model_mae_comparison.png`: across-seed model comparison.
