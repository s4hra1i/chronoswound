# Model card

## System

ChronosWound 0.2.0 contains three distinct components. They must not be conflated.

| Component | Data | Intended purpose | Valid claim |
|---|---|---|---|
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

Human leave-one-array-out accuracy was 83.3% with a 51.6–97.9% exact 95% interval. Injured-only accuracy was 88.9% with a 51.8–99.7% interval. These intervals demonstrate severe statistical uncertainty.

Ridge cross-severity rat MAE was 22.4 hours for mild-to-severe transfer and 31.7 hours for severe-to-mild transfer, against a 42.4-hour median baseline. Random-forest MAEs were 32.4 and 37.3 hours; its latter direction had negative R², demonstrating material large-error risk and the value of a simpler benchmark.

## Known performance gaps

There is no external human cohort, blinded assessment, calibration by demographic group, RNA-degradation challenge, batch-transfer study, tissue-site validation or comparison with forensic experts. Performance for post-mortem wounds is entirely unknown.

## Human oversight

No output should be interpreted without a qualified forensic pathologist and molecular scientist. At the present validation stage, such oversight does not make the tool suitable for casework; it only reduces the risk of misdescribing exploratory findings.
