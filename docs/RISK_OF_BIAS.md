# Risk-of-bias register

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
