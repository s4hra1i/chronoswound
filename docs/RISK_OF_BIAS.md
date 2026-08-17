# Risk-of-bias register

## Primary GSE178411 analysis

| Domain | Risk | Mitigation in this repository | Residual problem |
|---|---|---|---|
| Validation setting | High | Patient-grouped nested cross-validation and a prospectively tagged protocol | All evaluation remains internal resampling of one cohort |
| Sample size | High | All predictions, folds and patient assignments are retained | Only 49 samples from 39 patients |
| Repeated observations | Moderate | Every patient's samples remain in one fold | Only nine patients contribute repeated eligible samples; 31 are singletons |
| Predictor selection | Moderate | The 15-gene panel was fixed before GSE178411 outcome modelling | The biological panel was not independently selected or validated by another group |
| Outcome distribution | High | The exact target distribution and across-seed results are reported | The target contains only 16 distinct days and is concentrated around surgical scheduling |
| Late-range performance | High | Sample-level out-of-fold predictions and diagnostics are retained | Very few observations occur after day 14, producing large errors for some late wounds |
| Clinical confounding | High | Covariates-only and combined models are reported separately | Burn severity, treatment, infection, comorbidity and operative decision-making are incompletely measured |
| Storage and RNA quality | High | The limitation is explicitly reported | Samples span 2002–2018, but sample-level collection year and RNA integrity are unavailable |
| Normalisation and filtering | Moderate | CPM transformation is sample-local and preprocessing occurs within folds | The study authors' inherited count filtering cannot be reconstructed from the submitted matrix |
| Selective reporting | Low to moderate | The protocol, negative transfer results, full predictions and post-hoc labels are public | Repository provenance depends partly on the author's stated pre-analysis conduct |
| Transportability | High | Claims are restricted to surgically sampled burn wounds from 3–27 days | No independent site, platform, wound mechanism or post-mortem validation exists |
| Individual uncertainty | High | Population-level performance uncertainty is reported | The primary model does not provide calibrated prediction intervals for individual wounds |

## Exploratory GSE8056 analysis

| Domain | GSE8056 risk | Mitigation in this repository | Residual problem |
|---|---|---|---|
| Unit of analysis | High | Arrays are explicitly described as pooled samples | Patient-level error cannot be estimated |
| Sample size | High | Exact confidence intervals and all predictions are reported | Intervals remain very wide |
| Feature selection | Moderate | A locked biological panel is separated from exploratory genes | Initial panel was not formally preregistered |
| Control comparability | High | Injured-only results are reported separately | Normal tissue came from different procedures |
| Injury mechanism | High | Claims are restricted to thermal injury | Generalisation to blunt, sharp or post-mortem injury is unknown |
| Anatomy | High | Tissue-site variability is documented | Pool-level design prevents adjustment |
| Patient factors | High | Available pooled age and sex information is acknowledged | Individual covariates are unrecoverable |
| Batch/platform | Moderate | One platform and processed matrix are used consistently | No independent human platform validation |
| Outcome precision | High | Classification uses published time bins | Exact injury time cannot be inferred |
| Selective reporting | Moderate | Machine-readable full rankings and failures are retained | Only public deposited variables can be assessed |

## Cross-context dataset

GSE162565 improves the unit of analysis because its 33 arrays correspond to individual animals and it includes exact 1, 3, 24, 48 and 168-hour points under mild and severe injury. Its relevance to human skin remains indirect: it is rat skeletal muscle subjected to controlled contusion. It tests robustness across severity and biological context, not clinical transportability.
