# Changelog

## 0.3.1 — 2026-08-17

- Completed the planned 50-sample panel-only sensitivity analysis.
- Updated the Streamlit explorer to include the primary human analysis.
- Added late-range error diagnostics from stored out-of-fold predictions.
- Added GSE178411 to the dataset inventory and risk-of-bias register.
- Reconciled the results manuscript with the primary prospective analysis.
- Clarified that primary uncertainty estimates concern population-level performance rather
  than calibrated individual prediction intervals.

## 0.3.0 — 2026-08-17

- Published the GSE178411 analysis protocol before expression–outcome modelling.
- Added audited metadata for all 108 samples and a reproducible SOFT parser.
- Implemented the registered four-model repeated patient-grouped nested evaluation.
- Retained all folds, predictions, assignments, BCa uncertainty and 199 full-pipeline
  patient-block permutations.
- Reported the positive primary result with its fixed range, cohort and confounding limits.

## 0.2.1 — 2026-08-16

- Correct split-conformal calibration so the deployed model is not refitted after calibration.
- Add finite-sample quantile correction and report interval width alongside coverage.
- Add deterministic calibration and cross-context model tests.
- Use "interpretable" rather than overstating the scope of model explanation.

## 0.2.0 — 2026-08-15

- Added reproducible GSE8056 human burn-wound transcriptomics analysis.
- Added locked-panel leave-one-array-out evaluation, exact confidence intervals and permutation testing.
- Added GSE162565 cross-severity validation on individual rat muscle contusions.
- Added dataset card, study protocol, biomarker protocol and risk-of-bias register.
- Added an optional Streamlit research explorer.

## 0.1.0 — 2026-08-15

- Initial synthetic multimodal modelling pipeline with grouped validation and prediction intervals.
