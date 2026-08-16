"""Model training, grouped evaluation and uncertainty estimation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.compose import TransformedTargetRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, GroupShuffleSplit, cross_val_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import RobustScaler

from .data import GENES, validate_dataset
from .features import PhaseScoreTransformer

BASE_FEATURES = [*GENES, "neutrophil_pct", "fibroblast_pct", "temperature_c"]


@dataclass
class TrainingResult:
    model: Pipeline
    predictions: pd.DataFrame
    metrics: dict
    importance: pd.DataFrame
    interval_radius: float


def _candidate_models(seed: int) -> dict[str, Pipeline]:
    elastic = Pipeline([
        ("phase_scores", PhaseScoreTransformer()),
        ("impute", SimpleImputer(strategy="median")),
        ("scale", RobustScaler()),
        ("model", TransformedTargetRegressor(regressor=ElasticNet(alpha=0.08, l1_ratio=0.35, max_iter=10000), func=np.log1p, inverse_func=np.expm1)),
    ])
    forest = Pipeline([
        ("phase_scores", PhaseScoreTransformer()),
        ("impute", SimpleImputer(strategy="median")),
        ("model", RandomForestRegressor(n_estimators=350, min_samples_leaf=3, max_features=0.8, random_state=seed, n_jobs=-1)),
    ])
    return {"elastic_net": elastic, "random_forest": forest}


def train_and_evaluate(df: pd.DataFrame, seed: int = 42, coverage: float = 0.90) -> TrainingResult:
    validate_dataset(df)
    X = df[BASE_FEATURES]
    y = df["wound_age_hours"].to_numpy()
    groups = df["donor_id"].to_numpy()

    outer = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed)
    dev_idx, test_idx = next(outer.split(X, y, groups))
    X_dev, X_test = X.iloc[dev_idx], X.iloc[test_idx]
    y_dev, y_test = y[dev_idx], y[test_idx]
    groups_dev = groups[dev_idx]

    calibrator = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed + 1)
    train_rel, cal_rel = next(calibrator.split(X_dev, y_dev, groups_dev))
    cv = GroupKFold(n_splits=5)
    candidates = _candidate_models(seed)
    cv_mae = {}
    for name, estimator in candidates.items():
        scores = cross_val_score(estimator, X_dev.iloc[train_rel], y_dev[train_rel], groups=groups_dev[train_rel], cv=cv, scoring="neg_mean_absolute_error")
        cv_mae[name] = float(-scores.mean())

    winner = min(cv_mae, key=cv_mae.get)
    model = candidates[winner]
    model.fit(X_dev.iloc[train_rel], y_dev[train_rel])
    cal_pred = np.maximum(0, model.predict(X_dev.iloc[cal_rel]))
    residuals = np.abs(y_dev[cal_rel] - cal_pred)
    radius = float(np.quantile(residuals, coverage, method="higher"))

    model.fit(X_dev, y_dev)
    predicted = np.maximum(0, model.predict(X_test))
    lower = np.maximum(0, predicted - radius)
    upper = predicted + radius
    metrics = {
        "selected_model": winner,
        "cv_mae_hours": cv_mae,
        "test_mae_hours": float(mean_absolute_error(y_test, predicted)),
        "test_rmse_hours": float(mean_squared_error(y_test, predicted) ** 0.5),
        "test_r2": float(r2_score(y_test, predicted)),
        "target_interval_coverage": coverage,
        "observed_interval_coverage": float(np.mean((y_test >= lower) & (y_test <= upper))),
        "interval_radius_hours": radius,
        "n_development": int(len(dev_idx)),
        "n_test": int(len(test_idx)),
    }
    predictions = pd.DataFrame({"donor_id": groups[test_idx], "observed_hours": y_test, "predicted_hours": predicted, "lower_hours": lower, "upper_hours": upper})

    transformed = model[:-1].transform(X_test)
    names = model[:-1].get_feature_names_out() if hasattr(model[:-1], "get_feature_names_out") else [*BASE_FEATURES, "inflammation_score", "repair_score", "remodelling_score", "immune_balance"]
    pi = permutation_importance(model[-1], transformed, y_test, scoring="neg_mean_absolute_error", n_repeats=15, random_state=seed)
    importance = pd.DataFrame({"feature": names, "importance": pi.importances_mean}).sort_values("importance", ascending=False)
    return TrainingResult(model, predictions, metrics, importance, radius)
