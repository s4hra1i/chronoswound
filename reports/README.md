# Reproduced outputs

This directory contains the checked outputs used in the README and results manuscript.

- `gse8056/` contains the human exploratory analysis, locked-panel predictions, full temporal ranking and figures.
- `gse162565/` contains cross-severity predictions, model/baseline metrics and the comparison figure.
- `gse178411/` contains every output from the prospectively tagged primary human analysis,
  including all repeated out-of-fold predictions, folds, assignments, BCa comparison summary
  and the complete 199-permutation null distribution. It also contains the planned 50-sample
  panel-only sensitivity outputs and post-hoc sample-level error diagnostics.
- `synthetic/` contains the generated software demonstration, explicitly separated from all
  biological evidence. Its performance must not be reported as wound-age validation.

Every result can be regenerated from public GEO inputs using the commands in the repository README. Generated outputs are committed so reviewers can inspect the evidence without first downloading platform files. They are not additional independent datasets.
