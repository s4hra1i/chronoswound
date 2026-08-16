"""Persist model artefacts and diagnostic visualisations."""

from __future__ import annotations

import json
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import seaborn as sns

from .modelling import TrainingResult


def save_report(result: TrainingResult, output: str | Path) -> None:
    root = Path(output)
    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    (root / "metrics.json").write_text(json.dumps(result.metrics, indent=2) + "\n")
    result.predictions.to_csv(root / "predictions.csv", index=False)
    joblib.dump({"model": result.model, "interval_radius": result.interval_radius, "metrics": result.metrics}, root / "model.joblib")

    sns.set_theme(style="whitegrid", context="talk")
    p = result.predictions
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.errorbar(p["observed_hours"], p["predicted_hours"], yerr=[p["predicted_hours"] - p["lower_hours"], p["upper_hours"] - p["predicted_hours"]], fmt="o", alpha=.55, color="#7c3aed", ecolor="#c4b5fd")
    lim = max(p["observed_hours"].max(), p["upper_hours"].max())
    ax.plot([0, lim], [0, lim], "--", color="#334155", label="Ideal")
    ax.set(xlabel="Observed wound age (hours)", ylabel="Predicted wound age (hours)", title="Predictions with 90% intervals")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figures / "predicted_vs_observed.png", dpi=180)
    plt.close(fig)

    residual = p["predicted_hours"] - p["observed_hours"]
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.scatterplot(x=p["predicted_hours"], y=residual, ax=ax, color="#0f766e", alpha=.7)
    ax.axhline(0, ls="--", color="#334155")
    ax.set(xlabel="Predicted age (hours)", ylabel="Residual (hours)", title="Residual diagnostic")
    fig.tight_layout()
    fig.savefig(figures / "residuals.png", dpi=180)
    plt.close(fig)

    top = result.importance.head(12).sort_values("importance")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.barh(top["feature"], top["importance"], color="#f59e0b")
    ax.set(xlabel="Increase in MAE after permutation", title="Hold-out feature importance")
    fig.tight_layout()
    fig.savefig(figures / "feature_importance.png", dpi=180)
    plt.close(fig)
