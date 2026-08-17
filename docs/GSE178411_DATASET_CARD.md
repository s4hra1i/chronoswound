# GSE178411 dataset card

## Identity and source

- **Accession:** GSE178411
- **GEO title:** *Whole-transcriptome analysis illustrates evolving transcriptional human
  response to injury in acute wounds and scars*
- **Organism:** *Homo sapiens*
- **Assay:** bulk RNA sequencing on Illumina HiSeq instruments
- **Reference genome/annotation:** GRCh38 with the Rsubread inbuilt gene annotation
- **GEO record:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178411
- **Public date:** 3 June 2024
- **Repository retrieval/audit date:** 16 August 2026
- **Linked publication:** no peer-reviewed primary article was linked from the GEO series record
  when this card was prepared. The dataset is therefore cited directly rather than assigned an
  inferred publication.

## Original study

The submitted series contains 108 de-identified human skin samples: uninjured skin, acute burn
wounds and hypertrophic scars. Samples were collected from people undergoing operations at the
UW Medicine Regional Burn Center between 2002 and 2018. The GEO record states that tissue was
snap-frozen or placed in RNAlater before long-term storage at -80°C.

The processed matrix contains study-author-generated raw gene counts. The submitters used
Rsubread 1.32.4 and filtered genes using edgeR 3.24.3 `filterByExpr`. Raw reads are available
through the linked SRA records, but ChronosWound uses the submitted processed count matrix.

## ChronosWound cohort

The registered primary comparison retains samples that:

1. are labelled as wound tissue;
2. are labelled `Early Wound` or `Late wound`;
3. have a numeric day since injury; and
4. have recorded age for the four-model comparison.

This produces 49 samples from 39 patients spanning 3–27 days. The planned molecular-panel-only
sensitivity analysis includes `GSM5390619`, whose age is missing, and therefore contains 50
samples from 40 patients. Chronic wounds at 56, 97 and 1,291 days are outside the registered
estimand and are excluded.

The committed sample metadata are derived from GEO sample records by
`scripts/parse_gse178411_metadata.py`. The analysis downloads the submitted count matrix from
NCBI and verifies SHA-256 checksum
`19622a1b543d9b67481ca5bb13e35f73aafd28fe56f62dde3309e91c37ad0228`.

## Processing used here

ChronosWound selects the 15 locked Entrez Gene IDs, divides each count by that sample's complete
submitted library size, multiplies by one million and applies `log2(CPM + 0.5)`. This operation
is sample-local. Encoding, imputation, scaling, ridge fitting and hyperparameter selection are
performed within training data.

## Known limitations

- This is one centre's cohort and does not constitute external validation.
- Only nine eligible patients have repeated samples; 31 are singletons.
- The target is concentrated at surgical scheduling intervals and only seven primary samples
  occur after day 14.
- Burn severity, treatment, infection and other potential confounders are incompletely recorded.
- Sample-level collection year and RNA-integrity measurements are unavailable, despite collection
  spanning 2002–2018.
- Long-term storage and degradation may therefore be associated with apparent wound stage.
- The submitted matrix has already undergone study-author `filterByExpr` filtering, which cannot
  be reconstructed exactly from the processed matrix alone.
- GEO does not state a general-purpose licence for the deposited records. Users must follow NCBI
  terms and assess reuse requirements for their own application.

## Appropriate use

This dataset supports a narrowly scoped internal evaluation of temporal molecular signal in
surgically sampled subacute burn wounds. It must not be used to claim validation for wounds under
72 hours, other mechanisms, post-mortem tissue, individual prediction intervals, clinical use or
forensic casework.
