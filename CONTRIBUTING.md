# Contributing

Contributions are welcome when they improve reproducibility, validation or biological interpretation.

1. Open an issue describing the scientific or software change.
2. Create a focused branch and add tests for changed behaviour.
3. Run `pytest` and `ruff check src tests`.
4. Do not commit participant-level data, generated GEO downloads or fitted models containing sensitive information.
5. Report negative or null validation results with the same prominence as positive ones.

New biomarkers must be proposed with a citation and biological rationale before they are evaluated on a held-out dataset. Post-hoc additions belong in exploratory analyses and must not be described as confirmatory.
