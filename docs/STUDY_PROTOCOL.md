# Analysis protocol

## Primary question

Can a small, biologically fixed transcript panel distinguish published wound-age intervals without selecting genes on the held-out array?

## Secondary questions

1. Does performance persist after normal controls are removed?
2. Are errors at least confined to adjacent temporal bins?
3. Can the same broad marker panel estimate time when injury severity changes in an independent animal experiment?

## Outcomes

- Primary: leave-one-array-out balanced accuracy across four GSE8056 groups.
- Key sensitivity outcome: balanced accuracy across the three injured time bins.
- Uncertainty: exact 95% binomial interval and a 999-permutation null test.
- Cross-context outcomes: MAE, bootstrap 95% interval, RMSE, R² and proportion within 24 hours for fixed ridge and random-forest models when training on one contusion severity and testing the other; both are compared with a training-median baseline.

## Interpretation thresholds

No threshold converts these analyses into forensic validity. Results are considered hypothesis-generating even if statistically distinguishable from the permutation null. Failure on injured-only or cross-severity analysis outweighs strong normal-versus-injured separation.

## Reproducibility

All processed inputs are downloaded from NCBI GEO; probe mappings, analysis code, seeds, predictions and full gene-level results are retained. Generated files are not used as source inputs for subsequent model fitting.
