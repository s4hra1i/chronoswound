# Data

## GSE178411 metadata

`GSE178411_sample_metadata.csv` is a reproducible extraction of the public,
de-identified sample characteristics in the GEO SOFT family record. Regenerate it with:

```bash
python scripts/parse_gse178411_metadata.py \
  data/geo/GSE178411_family.soft.gz \
  data/GSE178411_sample_metadata.csv
```

Source records:

- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178411/soft/GSE178411_family.soft.gz`
- `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178411/suppl/GSE178411%5Fcounts.txt.gz`

The raw SOFT and count files are downloaded locally and are not committed. The metadata use
`wound_stage` for the first GEO `wound type` characteristic (for example, `Early Wound`) and
`sample_class` for the second (`wound`, `scar`, or `uninjured`).

The intended primary cohort is not silently encoded in the CSV. Its eligibility rules will be
timestamped separately before modelling.

## Synthetic fixture

Generate the synthetic software-test fixture with:

```bash
chronoswound generate --output data/synthetic_wounds.csv
```

The generated values are circular by design and test software behaviour only. They must not be
reported as evidence of biological or predictive performance.
