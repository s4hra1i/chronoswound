# Rapid literature and dataset search

## Scope

This is a reproducible rapid scoping search, not a systematic review. It was performed to define a biologically plausible panel and locate public data suitable for honest validation.

## Search date and sources

- Last searched: 15 August 2026
- Sources: NCBI PubMed/PMC, NCBI Gene Expression Omnibus and citation chaining from included primary studies

## Search concepts

Searches combined the following concepts:

```text
(wound age OR injury age OR post-injury interval)
AND (gene expression OR transcriptomic OR microarray OR RNA-seq)
AND (skin OR muscle OR wound)
```

GEO-specific searches included:

```text
"wound" AND "time course" AND "Homo sapiens"
"wound age estimation" AND "expression profiling"
"burn wound margin" AND "post burn day"
```

## Dataset inclusion criteria

- Transcript-level measurements from wounded tissue
- Known post-injury times or intervals
- At least two injured time points
- Processed public expression data and sufficient sample labels
- In vivo tissue preferred over isolated-cell scratch assays

## Exclusion criteria

- Treatment-only comparisons without injury timing
- Chronic ulcers without a known injury onset
- Cell migration assays described as “wound healing” but lacking actual tissue injury
- Studies without reusable expression matrices
- Datasets whose sample annotations could not distinguish time groups

## Included evidence

| Study | Relevance | Decision |
|---|---|---|
| Greco *et al.*, GSE8056 | Human burn-margin tissue; 0–3, 4–7 and >7-day pools | Primary human exploratory dataset |
| Li *et al.*, GSE162565 | Individual rat muscle contusions; two severities; 1–168 hours | Cross-severity and cross-context validation |
| Greco *et al.* primary paper | Biological interpretation of GSE8056 | Included as primary source |
| Li *et al.* primary paper | Design and limitations of GSE162565 | Included as primary source |

## Marker-selection approach

The panel was chosen by biological process rather than by searching validation results for the most accurate genes. It covers inflammatory signalling, leukocyte recruitment, matrix degradation, angiogenesis, matrix deposition and re-epithelialisation. Exact rules are recorded in `BIOMARKER_PROTOCOL.md`.

## Remaining gap

The search did not identify a suitable independent, unpooled human cohort containing both expression data and reliable post-injury timing. GSE162565 improves sample independence and time resolution but changes species, tissue and injury mechanism. A future formal systematic search should use at least two independent screeners, database-specific syntax, duplicate removal, a PRISMA flow diagram and risk-of-bias assessment at study level.
