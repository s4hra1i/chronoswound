# Dataset card: GSE8056

## Intended use

GSE8056 is used for an exploratory, reproducible demonstration of temporal gene-expression analysis in injured human skin. It is **not** used to train or validate the exact-hour ChronosWound estimator.

## Study design

- Organism: *Homo sapiens*
- Tissue: viable margins of thermally injured skin; uninjured human skin controls
- Assay: Affymetrix Human Genome U133 Plus 2.0 microarray (GPL570)
- Groups: 0–3 days, 4–7 days, >7 days after thermal injury, and normal skin
- Arrays: three per group, twelve total
- Pooling: each array combines equal masses of RNA from five tissue specimens
- Processing reported by depositors: RMA in GeneSpring 7.0

Source: Greco III JA *et al.* Gene expression analysis of injured human skin reveals dramatic early upregulation of extracellular matrix-related genes. *Journal of Investigative Dermatology* (2010). GEO accession: GSE8056; PMID: 19781859.

## Critical limitations

The array is the unit of analysis, not the individual patient. Pooling prevents donor-level uncertainty estimation and conceals within-group heterogeneity. The intervals are coarse, injury mechanism is restricted to burns, tissue sites vary, and the control tissue came from different surgical procedures. With only three arrays per group, p-values and FDR rankings are exploratory. No classifier or exact-time performance claim from this dataset would be credible without independent validation.

## Reproducibility

`chronoswound real-analysis` downloads the processed series matrix and GPL570 annotation directly from NCBI GEO, records sample labels, maps probes to gene symbols, log2-transforms values, collapses duplicate probes by the median, and generates PCA, marker-panel and temporal-signal outputs.
