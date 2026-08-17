# Model card

## System

ChronosWound 0.3.1 contains four distinct components. They must not be conflated.

| Component | Data | Intended purpose | Valid claim |
|---|---|---|---|
| Prospective human regressor | 49 GSE178411 wound samples from 39 patients | Test exact-day resolution from 3–27 days | Internally validated temporal signal in this cohort and range |
| Human fixed-panel classifier | 12 pooled GSE8056 arrays | Test broad temporal separability | Exploratory array-level time-bin classification |
| Rat cross-severity regressor | 30 individual GSE162565 animals | Test sensitivity to injury severity | Cross-context robustness evidence |
| Synthetic multimodal regressor | Generated observations | Demonstrate software and uncertainty workflow | Engineering demonstration only |

## Prohibited uses

- Dating a real injury
- Supporting criminal, civil, clinical or insurance decisions
- Reporting the synthetic 6.6-hour MAE as biological performance
- Treating pooled human arrays as individual patients
- Converting broad GSE8056 bins into exact hours

## Evaluation summary

The prospective GSE178411 locked-panel ridge achieved mean repeated patient-grouped MAE of 2.80
days, versus 4.26 days for the training-fold median. The 34.3% reduction met the pre-tagged
success rule: paired patient-clustered BCa improvement 1.46 days (95% CI 0.68–2.43) and
full-pipeline patient-block permutation *p*=0.005. Covariates alone achieved 4.58 days; combining
them with the panel achieved 3.03 days, worse than the panel alone.

The planned 50-sample panel-only sensitivity analysis achieved MAE of 3.00 days versus 4.41
days for its training-fold median. Post-hoc sample-level diagnostics exposed sparse late-range
performance: the sole 27-day sample had a mean out-of-fold prediction of approximately 13.8
days, and individual 19- and 22-day samples had errors above eight days.

Human leave-one-array-out accuracy was 83.3% with a 51.6–97.9% exact 95% interval. Injured-only accuracy was 88.9% with a 51.8–99.7% interval. These intervals demonstrate severe statistical uncertainty.

Ridge cross-severity rat MAE was 22.4 hours for mild-to-severe transfer and 31.7 hours for severe-to-mild transfer, against a 42.4-hour median baseline. Random-forest MAEs were 32.4 and 37.3 hours; its latter direction had negative R², demonstrating material large-error risk and the value of a simpler benchmark.

## Known performance gaps

There is no external human validation, blinded assessment, calibration by demographic group,
RNA-degradation challenge, batch-transfer study, independent tissue-site validation or
comparison with forensic experts. GSE178411 provides no evidence below three days, and
performance for post-mortem wounds is entirely unknown.

The primary GSE178411 uncertainty estimates apply to cohort-level performance comparisons. They
are not calibrated prediction intervals for individual wounds.

## Human oversight

No output should be interpreted without a qualified forensic pathologist and molecular scientist. At the present validation stage, such oversight does not make the tool suitable for casework; it only reduces the risk of misdescribing exploratory findings.
