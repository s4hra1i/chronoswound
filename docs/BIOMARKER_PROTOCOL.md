# Locked biomarker protocol

## Purpose

This protocol defines the marker panel before validation analyses are interpreted. Its purpose is to reduce outcome-driven gene selection and make failures visible.

## Inclusion logic

Markers were selected to cover complementary biological processes expected across wound repair:

| Process | Locked markers | Rationale |
|---|---|---|
| Acute inflammatory signalling | `IL6`, `TNF`, `CXCL8`, `CXCL2`, `PTGS2` | Cytokine and chemokine recruitment following tissue injury |
| Neutrophil/myeloid response | `MPO`, `CD68`, `CCL2` | Early granulocyte activity and later monocyte/macrophage recruitment |
| Matrix degradation | `MMP9` | Extracellular-matrix turnover during inflammation and repair |
| Angiogenesis and repair | `VEGFA`, `TGFB1` | Vascular response and transition towards repair |
| Matrix deposition | `COL1A1`, `FN1`, `SERPINE1` | Fibroplasia, matrix organisation and protease regulation |
| Epidermal response | `KRT14` | Basal keratinocyte activation and re-epithelialisation |

## Locked analysis rules

1. All available unambiguous probes for a marker are collapsed by median expression.
2. No gene is added because it ranks highly in the validation dataset.
3. Human time-bin testing uses nearest-centroid classification, standardised inside each training fold.
4. Every human array is held out once; feature selection is not repeated because the panel is fixed in advance.
5. Four-group performance and the harder injured-only three-bin performance are both reported.
6. Cross-context analysis uses the orthologous rat markers available on GPL17117. `CXCL8` is absent because rodents have no direct one-to-one orthologue; no post-hoc substitute is introduced.
7. Cross-severity validation trains on mild contusions and tests severe contusions, then reverses the direction.

## Exclusions

Histamine and serotonin are not gene-expression measurements and therefore are not treated as transcript markers. Candidate genes found through the exploratory ANOVA are reported separately and never passed into the locked validation model.

## Status

This protocol is retrospective relative to the original independent-project idea but prospective relative to the added validation outputs. A future confirmatory study should timestamp and register an updated protocol before accessing its test cohort.
