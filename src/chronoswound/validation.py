"""Small-sample validation designed to expose, rather than conceal, uncertainty."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import beta
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix
from sklearn.model_selection import LeaveOneOut
from sklearn.neighbors import NearestCentroid
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .real_data import FOCUS_GENES, GROUP_ORDER


def _binomial_interval(successes: int, total: int, alpha: float = 0.05) -> list[float]:
    lower = 0.0 if successes == 0 else float(beta.ppf(alpha / 2, successes, total - successes + 1))
    upper = 1.0 if successes == total else float(beta.ppf(1 - alpha / 2, successes + 1, total - successes))
    return [lower, upper]


def _loocv_predict(X: pd.DataFrame, y: np.ndarray) -> np.ndarray:
    predictions = np.empty(len(y), dtype=object)
    for train, test in LeaveOneOut().split(X):
        model = make_pipeline(StandardScaler(), NearestCentroid())
        model.fit(X.iloc[train], y[train])
        predictions[test[0]] = model.predict(X.iloc[test])[0]
    return predictions


def _permutation_p_value(X: pd.DataFrame, y: np.ndarray, observed: float, seed: int) -> float:
    rng = np.random.default_rng(seed)
    null = np.empty(999)
    for i in range(len(null)):
        shuffled = rng.permutation(y)
        null[i] = balanced_accuracy_score(shuffled, _loocv_predict(X, shuffled))
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


def evaluate_human_time_bins(
    genes: pd.DataFrame, output: str | Path, seed: int = 42
) -> dict:
    """Evaluate a fixed marker panel with honest leave-one-array-out testing."""
    root = Path(output)
    root.mkdir(parents=True, exist_ok=True)
    (root / "figures").mkdir(exist_ok=True)
    present = [g for g in FOCUS_GENES if g in genes.index]
    X = genes.loc[present].T
    y = np.asarray(["0–3 days"] * 3 + ["4–7 days"] * 3 + [">7 days"] * 3 + ["Normal"] * 3)
    pred = _loocv_predict(X, y)
    accuracy = float(accuracy_score(y, pred))
    balanced = float(balanced_accuracy_score(y, pred))

    injured = y != "Normal"
    injured_pred = _loocv_predict(X.loc[injured], y[injured])
    injured_accuracy = float(accuracy_score(y[injured], injured_pred))
    injured_balanced = float(balanced_accuracy_score(y[injured], injured_pred))
    ranks = {label: i for i, label in enumerate(GROUP_ORDER[1:])}
    adjacent = float(
        np.mean([abs(ranks[a] - ranks[b]) <= 1 for a, b in zip(y[injured], injured_pred)])
    )

    predictions = pd.DataFrame(
        {"sample": genes.columns, "observed": y, "predicted": pred, "correct": y == pred}
    )
    predictions.to_csv(root / "human_loocv_predictions.csv", index=False)
    summary = {
        "feature_strategy": "fixed literature-informed panel",
        "n_markers": len(present),
        "markers": present,
        "all_four_groups_accuracy": accuracy,
        "all_four_groups_accuracy_95ci": _binomial_interval(int((y == pred).sum()), len(y)),
        "all_four_groups_balanced_accuracy": balanced,
        "all_four_groups_permutation_p": _permutation_p_value(X, y, balanced, seed),
        "injured_only_three_bins_accuracy": injured_accuracy,
        "injured_only_accuracy_95ci": _binomial_interval(
            int((y[injured] == injured_pred).sum()), int(injured.sum())
        ),
        "injured_only_balanced_accuracy": injured_balanced,
        "injured_only_adjacent_bin_accuracy": adjacent,
        "interpretation": "exploratory pooled-array classification; not patient-level validation",
    }
    (root / "human_validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    matrix = confusion_matrix(y, pred, labels=GROUP_ORDER)
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(matrix, annot=True, fmt="d", cmap="Purples", cbar=False, ax=ax)
    ax.set_xticklabels(GROUP_ORDER, rotation=30, ha="right")
    ax.set_yticklabels(GROUP_ORDER, rotation=0)
    ax.set(xlabel="Predicted", ylabel="Observed", title="Fixed-panel leave-one-array-out validation")
    fig.tight_layout()
    fig.savefig(root / "figures" / "human_loocv_confusion.png", dpi=180)
    plt.close(fig)
    return summary
