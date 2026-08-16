import numpy as np
import pandas as pd

from chronoswound.cross_context import rat_metadata
from chronoswound.validation import _binomial_interval, _loocv_predict


def test_exact_interval_contains_observed_proportion():
    lower, upper = _binomial_interval(10, 12)
    assert lower < 10 / 12 < upper


def test_loocv_returns_one_prediction_per_observation():
    X = pd.DataFrame({"a": [0, 0.1, 3, 3.1, 6, 6.1], "b": [0, 0, 1, 1, 2, 2]})
    y = np.array(["early", "early", "middle", "middle", "late", "late"])
    prediction = _loocv_predict(X, y)
    assert len(prediction) == len(y)


def test_rat_metadata_matches_deposited_design():
    meta = rat_metadata(pd.Index([f"S{i}" for i in range(33)]))
    assert len(meta) == 33
    assert (meta.severity == "control").sum() == 3
    assert set(meta.loc[meta.severity != "control", "hours"]) == {1, 3, 24, 48, 168}
