# Prospective GSE178411 analysis protocol

**Protocol version:** 1.0  
**Protocol seed:** `20260816`  
**Planned tag:** `gse178411-protocol-v1.0`

## Prospective boundary

Before this protocol was tagged, only the structure of the public metadata, cohort
eligibility, target distribution, repeated-patient structure and presence of the locked
markers in the submitted count matrix were audited. No expression–outcome association or
predictive model using GSE178411 was examined. Results will be committed separately and will
identify this protocol tag.

The earlier marker panel was documented as fixed within the repository, but not independently
verifiable from the commit history; the GSE8056 result is therefore reported as exploratory.

## Question and estimand

The primary question is whether a locked 15-gene transcript panel provides practically useful
resolution of wound age in surgically sampled burn wounds collected 3–27 days after injury.
This is an early- versus late-subacute burn-wound problem, not general time-since-injury
estimation and not evidence for wounds sampled before 72 hours, other injury mechanisms or
forensic casework.

## Cohort fixed before modelling

The public GSE178411 metadata contain 108 samples. The primary eligible cohort requires:

1. `sample_class == "wound"`;
2. `wound_stage` equal to `Early Wound` or `Late wound`;
3. a numeric `days_since_injury`; and
4. a recorded age for the four-model comparison.

This produces 49 samples from 39 patients over 3–27 days. `GSM5390619` is excluded from the
four-model comparison because age is missing. A panel-only sensitivity analysis includes that
sample, giving 50 samples from 40 patients.

The three `Chronic wound` samples at 56, 97 and 1,291 days are excluded because they represent
a different biological and clinical quantity, not because of an outcome-derived statistical
outlier rule.

Only nine of 40 patients contribute repeated eligible samples. Patient grouping prevents
leakage for those patients but will have limited practical influence for the 31 singleton
patients.

## Locked predictors

The molecular panel is `IL6`, `TNF`, `CXCL8`, `CXCL2`, `PTGS2`, `MPO`, `CD68`, `CCL2`,
`MMP9`, `VEGFA`, `TGFB1`, `COL1A1`, `FN1`, `SERPINE1` and `KRT14`. All 15 corresponding
NCBI Gene IDs occur in the submitted matrix after the study authors' prior `filterByExpr`
processing. No marker will be added or removed after outcome modelling.

Covariates are age, sex, burn type and anatomical location. Unknown values remain explicit
categorical levels. Race is described but excluded from prediction because several levels have
very small counts and its use would be difficult to justify at this sample size.

## Prespecified models

All models use identical outer folds and the 49 complete observations:

1. training-fold median;
2. covariates-only ridge;
3. locked-panel ridge; and
4. covariates-plus-panel ridge.

The regression target is untransformed days. Ridge alpha is selected only within each outer
training set using inner patient-grouped cross-validation over the exact grid
`10 ** (-3 + 0.25 * k)` for `k = 0, ..., 24`, spanning 0.001–1,000.

## Preprocessing boundary

Library-size CPM and `log2(CPM + 0.5)` are sample-local and may be computed before splitting.
No TMM factor, cross-sample normalisation, outcome-informed filtering or global gene-wise
centring will be fitted on the complete cohort. Categorical encoding and feature
standardisation are fitted inside each training fold and applied unchanged to its test fold.
The study authors' inherited count filtering is disclosed and cannot be reconstructed from the
submitted filtered matrix.

## Repeated grouped evaluation

The outer evaluation uses five-fold `GroupKFold(shuffle=True)` repeated for seeds `20260816`
through `20260835`. This requires scikit-learn 1.6 or later. Each patient remains wholly within
one fold. Inner grouped tuning uses a random state deterministically derived from the outer seed
and fold number.

The primary point estimate is the mean of the 20 repetition-specific pooled out-of-fold MAEs.
The report will include all 20 pooled MAEs, their median, IQR and range; all 100 individual-fold
MAEs; every out-of-fold prediction; and every patient-to-fold assignment.

All outcome units are days. The descriptive eight-day global-median MAE is 4.08 days on the
49-sample complete-case cohort, but the formal comparator is the median fitted within each
outer training fold.

## Prespecified success rule

Practically useful standalone molecular resolution is declared only if all three conditions
hold:

1. mean repeated-CV panel MAE is at least 20% lower than the corresponding training-median MAE;
2. the 95% patient-clustered BCa interval for `baseline MAE − panel MAE` lies wholly above zero;
3. the patient-block permutation test gives `p < 0.05`.

This threshold is deliberately conservative. Failure to clear it does not establish absence of
molecular signal; it establishes only that practically useful temporal resolution was not
demonstrated at this sample size.

Incremental molecular value beyond recorded covariates additionally requires the combined
model to beat the covariates-only model with a patient-clustered BCa interval for the paired MAE
difference wholly above zero. Smaller improvements may be reported but not described as
established incremental value.

## Uncertainty and null procedure

The BCa interval uses 10,000 bootstrap replicates of the paired, fixed out-of-fold absolute
errors, resampling patients and retaining all observations belonging to each selected patient.
Acceleration uses leave-one-patient-out jackknife statistics. This interval represents patient
sampling uncertainty conditional on the fitted repeated-CV procedure; it does not incorporate
new model fitting. A replicate with fewer than two unique patients or an undefined statistic is
rejected and redrawn.

The permutation test instead measures whether the complete modelling procedure manufactures
apparent skill under a null association. Day-value blocks are permuted between patients within
equal sample-count strata, preserving stable within-patient row order. The 31 singleton
patients therefore permute like individual samples, the eight two-sample patients permute as
blocks, and the sole three-sample patient cannot move. The stratification consequently has only
a small practical effect in this cohort, but avoids breaking the repeated-patient structure that
does exist.

For each of 199 permutations, the full 20-repetition nested grouped pipeline is rerun. The
one-sided Monte Carlo p-value is `(1 + number of null improvements at least as large as the
observed improvement) / 200`, giving minimum resolution 0.005. The
bootstrap and permutation address different uncertainty: the former patient sampling
conditional on fitted predictions, the latter the full procedure under a null outcome mapping.

## Secondary classification

Early/Late classification is descriptive only. The labels are a deterministic threshold of the
target: Early covers 3–7 days and Late 8–27 days. No success criterion attaches to this analysis,
and classification performance cannot replace or rescue a failed regression result.

## Known interpretive constraints

The target occurs on only 16 distinct days and is concentrated around surgical scheduling:
72% of eligible samples fall within 3–10 days. Burn severity, anatomical location and clinical
decision-making may influence operation timing, so covariate predictiveness may represent
scheduling rather than biological ageing.

Samples were collected between 2002 and 2018, but sample-level collection year and RNA-integrity
measurements are unavailable. RNA degradation during long-term storage could therefore correlate
with collection era and apparent wound stage, producing temporal signal that cannot be tested or
adjusted here.

## Reproducibility record

The results JSON will record the protocol tag, commit hash, Python, NumPy, pandas and
scikit-learn versions, all seeds, alpha grid, fold assignments, eligibility exclusions, model
metrics, bootstrap settings and permutation settings. Results—positive or negative—will be
committed after this protocol tag in a separate commit that references it.
