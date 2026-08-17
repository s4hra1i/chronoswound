# ChronosWound: prospective internal evaluation of a locked transcript panel for subacute burn-wound age

## Abstract

ChronosWound evaluates whether a biologically fixed transcript panel contains temporal
information after injury while making small-sample failure modes explicit. The primary analysis
used 49 surgically sampled human burn wounds from 39 patients in GSE178411, spanning 3–27 days
after injury. Its protocol, cohort rules, predictors, models and success criterion were tagged
before expression–outcome modelling. A locked 15-gene ridge model achieved mean repeated
patient-grouped out-of-fold mean absolute error (MAE) of 2.80 days, compared with 4.26 days for
the training-fold median. The paired improvement was 1.46 days (patient-clustered BCa 95% CI
0.68–2.43), and none of 199 full-pipeline patient-block permutations matched the observed result
(one-sided Monte Carlo *p*=0.005). A planned panel-only sensitivity analysis including the
sample excluded for missing age contained 50 samples from 40 patients and produced panel MAE
3.00 days versus 4.41 days for the training-fold median. Earlier pooled human arrays and rat
cross-severity transfer are retained as exploratory evidence. The primary result demonstrates
internal temporal resolution in one cohort, not external, clinical or forensic validation.

## Introduction

Wound repair is not a single molecular clock. Inflammation, cell recruitment, angiogenesis,
matrix turnover and re-epithelialisation overlap, while anatomy, severity, treatment, infection,
health and sample handling can change their timing. An apparently accurate model may therefore
learn patient identity, batch structure or clinical scheduling rather than elapsed biological
time.

ChronosWound was designed around that problem. It separates a prospectively specified primary
human analysis from earlier exploratory human and animal work, keeps repeated observations from
the same patient together, fits preprocessing within training data and retains negative and
asymmetric transfer results.

## Methods

### Primary human cohort

GSE178411 contains 108 human skin RNA-seq samples. The registered complete-case comparison
included wound samples labelled Early or Late wound, with a numeric day since injury and
recorded age. This produced 49 samples from 39 patients over 3–27 days. The target occurred on
16 distinct days and was concentrated around surgical scheduling intervals. Only nine patients
contributed repeated eligible samples; the remaining 31 were singletons.

The locked panel comprised `IL6`, `TNF`, `CXCL8`, `CXCL2`, `PTGS2`, `MPO`, `CD68`, `CCL2`,
`MMP9`, `VEGFA`, `TGFB1`, `COL1A1`, `FN1`, `SERPINE1` and `KRT14`. Counts were converted using
sample-local library-size CPM followed by `log2(CPM + 0.5)`. No cohort-wide normalisation,
outcome-informed filtering or global gene-wise centring was fitted before splitting.

Four models used identical outer folds: a training-fold median, covariates-only ridge, locked
panel ridge and covariates-plus-panel ridge. Five-fold patient-grouped evaluation was repeated
over 20 fixed seeds. Ridge alpha selection, encoding and scaling occurred inside each outer
training set using inner grouped cross-validation.

The registered success rule required all of: at least 20% lower mean repeated-CV MAE than the
training-median baseline; a patient-clustered BCa 95% interval for the paired improvement wholly
above zero; and a one-sided full-pipeline patient-block permutation *p* below 0.05.

### Planned sensitivity analysis

The protocol also specified a panel-only sensitivity analysis including one wound excluded from
the four-model comparison because age was missing. This analysis used 50 samples from 40
patients, the same locked panel, seeds and grouped nested procedure, and compared the panel with
the training-fold median. No separate success criterion was specified.

### Exploratory tracks

GSE8056 contains twelve pooled human burn-margin microarrays in three injury-time bins and a
normal-skin group. Its array, not an individual patient, is the unit of analysis. GSE162565
contains 30 individually sampled wounded rats under mild and severe muscle contusion at five
exact times. These tracks test broad separability and sensitivity to biological context; neither
is external validation of the primary human estimator.

## Results

### Primary analysis

| Model | Mean pooled MAE | Range across 20 seeds |
|---|---:|---:|
| Training-fold median | 4.26 d | 4.08–4.63 d |
| Covariates only | 4.58 d | 4.30–5.05 d |
| Locked panel | **2.80 d** | 2.54–3.30 d |
| Covariates + panel | 3.03 d | 2.60–3.46 d |

The locked panel reduced MAE by 34.3%. Its paired improvement over the training-fold median was
1.46 days (patient-clustered BCa 95% CI 0.68–2.43), with a full-pipeline patient-block
permutation *p*=0.005. All three registered success conditions were met.

The preregistered covariates-only versus combined contrast favoured the combined model, but the
covariates-only model underperformed the median baseline. A more informative post-hoc comparison
showed that adding the recorded covariates worsened panel MAE from 2.80 to 3.03 days. These data
do not show that the available clinical covariates improve the panel.

The planned 50-sample panel-only sensitivity produced MAE of 3.00 days versus 4.41 days for the
training-fold median, an improvement of 1.41 days. It did not alter the registered 49-sample
primary estimand or result.

### Error distribution

Performance was not uniform across the observed range. Averaging each sample's panel prediction
across the 20 registered repetitions exposed regression towards the cohort centre. The sole
27-day sample had a mean prediction of approximately 13.8 days; individual samples at 19 and
22 days also had errors above eight days. Very few observations occurred after day 14, so the
overall MAE must not be read as comparable resolution throughout the full 3–27-day range.

### Exploratory evidence

The GSE8056 locked-panel classifier achieved 83.3% four-group array-level accuracy (exact 95% CI
51.6–97.9%) and 88.9% across injured-only bins (95% CI 51.8–99.7%). These intervals remain wide,
and panel ordering before that analysis cannot be independently verified from commit history.

Rat ridge transfer was asymmetric. Training on mild and testing severe contusions produced
22.4-hour MAE and R² 0.610; reversing the direction produced 31.7-hour MAE and R² 0.018. A random
forest had negative R² in the failed direction. This is evidence of context sensitivity, not
human transportability.

## Discussion

The registered primary result supports temporal information in the locked panel for surgically
sampled burn wounds in this cohort's 3–27-day range. It is stronger than the earlier pooled-array
evidence because it uses individual samples, exact days and patient identifiers. It nevertheless
remains internal resampling of one small cohort.

The principal threats are target concentration around operation schedules; sparse late-range
coverage; incomplete burn severity, treatment and infection information; and unavailable
sample-level collection year and RNA integrity despite collection spanning 2002–2018. Storage
or degradation structure could therefore masquerade as temporal biology. The inherited
study-author count filtering also cannot be reconstructed from the submitted filtered matrix.

The reported intervals quantify uncertainty in cohort-level performance comparisons conditional
on the retained repeated-CV predictions. The primary model does not provide calibrated
prediction intervals for individual wound-age estimates. The synthetic demonstration includes
such intervals only to test software behaviour and supplies no biological evidence.

## Conclusion

ChronosWound demonstrates useful internal temporal resolution from a locked transcript panel in
one subacute human burn cohort, while exposing large late-range errors and failed cross-context
transfer. It is a reproducible hypothesis-generating research project, not a validated wound-age
estimator, and must not be used in clinical, legal or forensic casework.

## References

1. National Center for Biotechnology Information. Gene Expression Omnibus series GSE178411: *Whole-transcriptome analysis illustrates evolving transcriptional human response to injury in acute wounds and scars*. Public 3 June 2024. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178411. No linked peer-reviewed primary article was identified on the GEO record as of 17 August 2026.
2. Greco JA III, Pollins AC, Boone BE, Levy SE, Nanney LB. A microarray analysis of temporal gene expression profiles in thermally injured human skin. *Burns*. 2010;36(2):192–204. PMID: 19781859. GEO: GSE8056.
3. Li N, Li C, Li D, *et al.* Identifying biomarkers for evaluating wound extent and age in the contused muscle of rats using microarray analysis: a pilot study. *PeerJ*. 2021;9:e12709. GEO: GSE162565.
