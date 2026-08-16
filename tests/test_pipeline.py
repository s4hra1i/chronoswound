import pandas as pd
import pytest

from chronoswound.data import generate_synthetic_cohort, validate_dataset
from chronoswound.features import PhaseScoreTransformer
from chronoswound.modelling import train_and_evaluate
from chronoswound.real_data import _benjamini_hochberg


def test_generator_is_reproducible():
    first = generate_synthetic_cohort(80, seed=7)
    second = generate_synthetic_cohort(80, seed=7)
    pd.testing.assert_frame_equal(first, second)


def test_phase_scores_are_added():
    df = generate_synthetic_cohort(80)
    transformed = PhaseScoreTransformer().fit_transform(df)
    assert {"inflammation_score", "repair_score", "remodelling_score", "immune_balance"} <= set(transformed)


def test_invalid_schema_is_rejected():
    df = generate_synthetic_cohort(80).drop(columns="IL6")
    with pytest.raises(ValueError, match="IL6"):
        validate_dataset(df)


def test_end_to_end_training():
    result = train_and_evaluate(generate_synthetic_cohort(140), seed=3)
    assert result.metrics["test_mae_hours"] > 0
    assert 0 <= result.metrics["observed_interval_coverage"] <= 1
    assert len(result.predictions) > 0


def test_fdr_adjustment_is_bounded_and_monotonic():
    adjusted = _benjamini_hochberg(pd.Series([0.001, 0.01, 0.2, 0.9]))
    assert adjusted.between(0, 1).all()
    assert adjusted.is_monotonic_increasing
