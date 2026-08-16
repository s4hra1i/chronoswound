# ChronosWound: an uncertainty-aware, cross-context analysis of temporal transcription after injury

## Abstract

Estimating the interval between injury and examination remains difficult because wound repair is dynamic, heterogeneous and affected by injury severity and patient characteristics. This project evaluated whether a small, biologically fixed transcript panel contains temporal information across two public microarray studies. GSE8056 comprised twelve human arrays representing pooled burn-margin samples collected 0–3, 4–7 or more than 7 days after injury and normal skin. A locked 15-gene panel spanning inflammatory recruitment, myeloid activity, angiogenesis, matrix turnover and re-epithelialisation was assessed without data-driven feature selection. Leave-one-array-out nearest-centroid classification achieved 83.3% accuracy across four groups (95% exact CI 51.6–97.9%; permutation *p*=0.001). After controls were removed, accuracy across the three injured intervals was 88.9% (95% CI 51.8–99.7%), although uncertainty remained substantial. Complementary analysis used GSE162565, comprising individual rat muscle contusions at 1, 3, 24, 48 and 168 hours. Training on mild injuries and testing severe injuries produced a mean absolute error (MAE) of 22.4 hours with ridge regression, compared with 32.4 hours for a random forest and 42.4 hours for the training-median baseline. Reversing the direction yielded MAEs of 31.7, 37.3 and 42.4 hours, respectively. These results support temporal signal in a compact injury-response panel but do not validate patient-level forensic estimation.

## Introduction

Wound age is not encoded by a single molecular clock. Haemostasis, inflammation, proliferation and remodelling overlap, while anatomical site, injury mechanism, severity, infection, medication and sampling conditions can alter their timing. Molecular approaches are nevertheless attractive because coordinated transcriptional changes may provide information beyond morphology alone.

Greco and colleagues profiled viable human burn margins and reported extensive expression changes, including prolonged inflammatory activity and extracellular-matrix perturbation. Their deposited GSE8056 dataset provides rare human temporal wound data but pools five specimens per array. Li and colleagues subsequently profiled individual rat muscle contusions at five exact post-injury times under two severities in GSE162565, explicitly examining biomarkers for wound extent and age. Together, these studies allow two distinct questions: whether a compact panel separates broad human burn intervals, and whether temporal prediction transfers when injury severity changes in a controlled model.

The project deliberately separates exploratory discovery, locked-panel evaluation and synthetic software demonstration. The primary aim was not to maximise accuracy, but to test whether interpretable temporal information survives validation procedures designed to expose small-sample uncertainty.

## Methods

### Human dataset and preprocessing

Processed RMA expression values and GPL570 annotations were downloaded directly from NCBI GEO. Probe intensities were transformed as log₂(*x*+1), mapped to unambiguous gene symbols and collapsed by median expression when multiple probes represented the same gene. The dataset contained three pooled arrays per group: 0–3 days, 4–7 days, >7 days and normal skin.

PCA used the 1,000 most variable genes after per-gene standardisation. Exploratory temporal differences were assessed by one-way ANOVA across the four groups, followed by Benjamini–Hochberg false-discovery-rate correction. These rankings were never used as validation features.

### Locked marker panel

The fixed panel represented acute cytokine signalling (`IL6`, `TNF`, `CXCL8`, `CXCL2`, `PTGS2`), neutrophil and myeloid activity (`MPO`, `CD68`, `CCL2`), matrix degradation (`MMP9`), angiogenesis and repair (`VEGFA`, `TGFB1`), matrix deposition (`COL1A1`, `FN1`, `SERPINE1`) and epidermal response (`KRT14`). Fourteen were available after unambiguous GPL570 mapping.

Nearest-centroid classification was assessed by leaving out each array once. Standardisation was fitted inside each training fold. Accuracy, balanced accuracy, an exact binomial 95% interval and a 999-permutation *p*-value were reported. A prespecified sensitivity analysis removed normal controls and classified only the three injury intervals.

### Cross-context analysis

GSE162565 contained three controls and 30 individually sampled wounded rats: three biological replicates at 1, 3, 24, 48 and 168 hours under mild and severe muscle contusion. Fourteen panel genes mapped unambiguously to GPL17117; `CXCL8` was omitted because rodents lack a direct one-to-one orthologue. Fixed ridge and random-forest regressors were trained on one severity and evaluated unchanged on the other. MAE, bootstrap 95% intervals, RMSE, R² and performance against a training-median baseline were retained.

## Results

### Human exploratory analysis

After annotation and probe collapsing, 20,848 genes across twelve arrays remained. PCA clearly separated normal skin from injured samples, while injured time groups showed partial rather than complete temporal ordering. The exploratory ANOVA identified 636 genes at FDR <0.05. This count should be interpreted cautiously because the arrays are pooled and the groups differ in more than elapsed time.

The locked panel correctly classified ten of twelve arrays (83.3%; exact 95% CI 51.6–97.9%). All three normal arrays, all three 0–3-day arrays and two of three 4–7-day and >7-day arrays were correctly assigned. Balanced accuracy was 83.3%, exceeding the permutation null (*p*=0.001).

Normal tissue did not solely drive the result. When controls were excluded and the model was refitted within each leave-one-out iteration, eight of nine wounded arrays were correctly classified (88.9%; exact 95% CI 51.8–99.7%). All injured-only predictions fell within the correct or an adjacent temporal interval. The very wide intervals remain more important than the point estimates.

### Cross-severity results

Training on mild contusions and testing severe contusions produced a ridge MAE of 22.4 hours (bootstrap 95% CI 8.6–41.1), RMSE of 38.7 hours and R² of 0.61. The random-forest MAE was 32.4 hours (95% CI 13.2–56.1), while the training-median baseline MAE was 42.4 hours. Thus the panel added temporal information in this direction, and the simpler regularised model transferred better.

Training on severe contusions and testing mild contusions produced a ridge MAE of 31.7 hours (bootstrap 95% CI 6.5–59.7), RMSE of 61.4 hours and R² of 0.02. The random forest reduced MAE from the 42.4-hour baseline to 37.3 hours but had negative R² (−0.24). This asymmetric result indicates that severity changes distort temporal profiles and that improved average absolute error can coexist with damaging large errors. Neither model can be described as severity-invariant.

## Discussion

The strongest finding is not a headline accuracy but the persistence of some temporal information under deliberately constrained analyses. Human classification remained high after controls were removed, reducing—but not eliminating—the concern that the model merely distinguished injured from uninjured tissue. Cross-severity rat testing also improved MAE in both directions, although the negative R² in one direction demonstrates poor transportability for some observations.

Biologically, the panel covers plausible transitions from cytokine and myeloid activation towards angiogenesis and matrix deposition. However, the analysis does not establish that any marker is specific to elapsed time. Changes may reflect burn depth, anatomical composition, systemic inflammation, age, sex or other pooled characteristics. The human array is the unit of analysis; treating its five component specimens as independent would be pseudoreplication.

The study has four central limitations. First, GSE8056 contains only twelve pooled arrays, yielding extremely wide confidence intervals. Second, its time labels are broad bins rather than exact outcomes. Third, normal tissue originated from different surgical settings. Fourth, the cross-context study changes species, tissue and injury mechanism simultaneously. It is valuable as a stress test, not external human validation.

A confirmatory study should use unpooled human samples, exact sampling times, predefined exclusion rules, RNA-quality and batch measures, and individual metadata covering anatomy, severity, infection, medication and comorbidity. The locked panel and model should then be frozen before testing at an independent site. Blinded comparison with forensic-pathologist estimates and explicit reporting of inconclusive results would be essential.

## Conclusion

ChronosWound shows that an interpretable injury-response panel contains temporal structure in pooled human burn data and retains limited signal across injury severity in an animal model. It also shows why apparent accuracy is insufficient: uncertainty is wide and cross-severity performance is asymmetric. The repository should therefore be understood as a reproducible hypothesis-generating framework, not a forensic instrument.

## References

1. Greco JA III, Pollins AC, Boone BE, Levy SE, Nanney LB. A microarray analysis of temporal gene expression profiles in thermally injured human skin. *Burns*. 2010;36(2):192–204. PMID: 19781859. GEO: GSE8056.
2. Li N, Li C, Li D, *et al.* Identifying biomarkers for evaluating wound extent and age in the contused muscle of rats using microarray analysis: a pilot study. *PeerJ*. 2021;9:e12709. GEO: GSE162565.
3. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society, Series B*. 1995;57(1):289–300.
