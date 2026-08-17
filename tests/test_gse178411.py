
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from chronoswound.gse178411 import (
    COUNTS_SHA256,
    PANEL,
    AnalysisConfig,
    _permuted_outcome,
    bca_cluster_interval,
    download_counts,
    repeated_nested_cv,
)


def _cohort() -> pd.DataFrame:
    rng = np.random.default_rng(7)
    rows = []
    for patient in range(15):
        day = 3 + patient
        row = {
            "geo_accession": f"GSM{patient}",
            "patient_id": patient,
            "days_since_injury": float(day),
            "age": 20 + patient,
            "sex": "female" if patient % 2 else "male",
            "burn_type": "flame" if patient % 3 else "contact",
            "location": "arm" if patient % 2 else "leg",
        }
        row.update({gene: day + rng.normal(0, 0.2) for gene in PANEL})
        rows.append(row)
    return pd.DataFrame(rows)


def test_repeated_nested_cv_keeps_patients_in_one_fold() -> None:
    cohort = _cohort()
    config = AnalysisConfig(
        seeds=(11, 12), outer_splits=3, inner_splits=3,
        bootstrap_replicates=50, permutation_replicates=1
    )
    predictions, folds, assignments = repeated_nested_cv(cohort, config)
    assert len(predictions) == len(cohort) * 4 * 2
    assert len(folds) == 3 * 4 * 2
    assert assignments.groupby(["repetition", "patient_id"])["fold"].nunique().max() == 1


def test_repeated_nested_cv_can_run_panel_only() -> None:
    cohort = _cohort()
    config = AnalysisConfig(
        seeds=(11,),
        outer_splits=3,
        inner_splits=3,
        bootstrap_replicates=50,
        permutation_replicates=1,
    )
    predictions, folds, _ = repeated_nested_cv(
        cohort, config, models=("panel",)
    )
    assert set(predictions["model"]) == {"median", "panel"}
    assert set(folds["model"]) == {"median", "panel"}


def test_cluster_bca_and_block_permutation() -> None:
    paired = pd.DataFrame(
        {
            "patient_id": np.repeat(np.arange(10), 2),
            "difference": np.linspace(0.1, 2.0, 20),
        }
    )
    interval = bca_cluster_interval(paired, replicates=200, seed=3)
    assert interval["lower_95_days"] < interval["estimate_days"] < interval["upper_95_days"]

    cohort = pd.DataFrame(
        {
            "patient_id": [1, 1, 2, 2, 3],
            "days_since_injury": [3.0, 5.0, 8.0, 10.0, 12.0],
        }
    )
    permuted = _permuted_outcome(cohort, np.random.default_rng(4))
    assert set(permuted[:4]) == {3.0, 5.0, 8.0, 10.0}
    assert permuted[4] == 12.0


def test_download_counts_rejects_unpinned_content(tmp_path) -> None:
    counts = tmp_path / "GSE178411_counts.txt.gz"
    counts.write_bytes(b"not the pinned count matrix")
    with pytest.raises(ValueError, match="checksum mismatch"):
        download_counts(counts)
    assert len(COUNTS_SHA256) == 64


def test_committed_primary_outputs_match_headline_metrics() -> None:
    root = Path(__file__).parents[1]
    predictions = pd.read_csv(root / "reports/gse178411/oof_predictions.csv")
    means = (
        predictions.groupby(["repetition", "model"])["absolute_error_days"]
        .mean()
        .groupby("model")
        .mean()
    )
    assert len(predictions) == 3920
    assert not predictions.duplicated(
        ["repetition", "model", "geo_accession"]
    ).any()
    assert means["median"] == pytest.approx(4.259183673469388)
    assert means["covariates"] == pytest.approx(4.578210842350088)
    assert means["panel"] == pytest.approx(2.799980146322827)
    assert means["combined"] == pytest.approx(3.0264050107599205)

    sensitivity = pd.read_csv(
        root / "reports/gse178411/sensitivity_repetition_metrics.csv"
    )
    sensitivity_means = sensitivity.groupby("model")["pooled_mae_days"].mean()
    assert sensitivity_means["median"] == pytest.approx(4.407)
    assert sensitivity_means["panel"] == pytest.approx(3.001254771632365)
