"""Cross-context validation using the rat muscle-contusion dataset GSE162565."""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.compose import TransformedTargetRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .real_data import parse_series_matrix

RAT_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162565/"
    "matrix/GSE162565_series_matrix.txt.gz"
)


def rat_metadata(columns: pd.Index) -> pd.DataFrame:
    severity = ["control"] * 3 + ["mild"] * 15 + ["severe"] * 15
    hours = [0] * 3 + [1] * 3 + [3] * 3 + [24] * 3 + [48] * 3 + [168] * 3
    hours += [1] * 3 + [3] * 3 + [24] * 3 + [48] * 3 + [168] * 3
    replicate = [1, 2, 3] * 11
    return pd.DataFrame(
        {"sample": columns, "severity": severity, "hours": hours, "replicate": replicate}
    )


def download_rat_matrix(data_dir: str | Path) -> Path:
    target = Path(data_dir) / "GSE162565_series_matrix.txt.gz"
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists():
        request = urllib.request.Request(RAT_MATRIX_URL, headers={"User-Agent": "ChronosWound/0.1"})
        with urllib.request.urlopen(request, timeout=90) as response, target.open("wb") as out:
            out.write(response.read())
    return target


def _metrics(y: np.ndarray, prediction: np.ndarray) -> dict:
    rng = np.random.default_rng(2026)
    bootstrap_mae = []
    for _ in range(2000):
        sample = rng.integers(0, len(y), len(y))
        bootstrap_mae.append(mean_absolute_error(y[sample], prediction[sample]))
    return {
        "mae_hours": float(mean_absolute_error(y, prediction)),
        "mae_bootstrap_95ci": [
            float(np.quantile(bootstrap_mae, 0.025)),
            float(np.quantile(bootstrap_mae, 0.975)),
        ],
        "rmse_hours": float(mean_squared_error(y, prediction) ** 0.5),
        "r2": float(r2_score(y, prediction)),
        "within_24_hours": float(np.mean(np.abs(y - prediction) <= 24)),
    }


def evaluate_cross_severity(
    matrix_path: str | Path,
    probe_map_path: str | Path,
    output: str | Path,
    seed: int = 42,
) -> dict:
    """Train on one injury severity and test unchanged on the other."""
    root = Path(output)
    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    probes = parse_series_matrix(matrix_path)
    mapping = pd.read_csv(probe_map_path).set_index("probe_id")["gene_symbol"]
    mapping.index = mapping.index.astype(str)
    probes.index = probes.index.astype(str)
    genes = probes.join(mapping.rename("gene"), how="inner")
    values = genes.drop(columns="gene")
    values["gene"] = genes["gene"]
    expression = values.groupby("gene").median().T
    metadata = rat_metadata(expression.index).set_index("sample")
    wounded = metadata.severity != "control"
    expression = expression.loc[wounded]
    metadata = metadata.loc[wounded]

    rows = []
    direction_metrics = {}
    for train_severity, test_severity in [("mild", "severe"), ("severe", "mild")]:
        train = metadata.severity == train_severity
        test = metadata.severity == test_severity
        model = RandomForestRegressor(
            n_estimators=1000, min_samples_leaf=2, max_features=0.7,
            random_state=seed, n_jobs=-1,
        )
        model.fit(expression.loc[train], np.log1p(metadata.loc[train, "hours"]))
        prediction = np.maximum(0, np.expm1(model.predict(expression.loc[test])))
        truth = metadata.loc[test, "hours"].to_numpy()
        ridge = make_pipeline(
            StandardScaler(),
            TransformedTargetRegressor(
                regressor=Ridge(alpha=10.0), func=np.log1p, inverse_func=np.expm1
            ),
        )
        ridge.fit(expression.loc[train], metadata.loc[train, "hours"])
        ridge_prediction = np.maximum(0, ridge.predict(expression.loc[test]))
        baseline = np.repeat(metadata.loc[train, "hours"].median(), len(truth))
        key = f"train_{train_severity}_test_{test_severity}"
        direction_metrics[key] = {
            "random_forest": _metrics(truth, prediction),
            "ridge": _metrics(truth, ridge_prediction),
            "training_median_baseline": _metrics(truth, baseline),
        }
        for model_name, model_prediction in [
            ("Random forest", prediction), ("Ridge", ridge_prediction)
        ]:
            for sample, observed, predicted_value in zip(
                metadata.index[test], truth, model_prediction
            ):
                rows.append(
                    {"sample": sample, "model": model_name,
                     "train_severity": train_severity,
                     "test_severity": test_severity, "observed_hours": observed,
                     "predicted_hours": predicted_value}
                )
    predictions = pd.DataFrame(rows)
    predictions.to_csv(root / "cross_severity_predictions.csv", index=False)
    summary = {
        "accession": "GSE162565",
        "organism": "Rattus norvegicus",
        "tissue": "skeletal muscle",
        "n_individual_wounded_animals": int(len(metadata)),
        "n_fixed_markers": int(expression.shape[1]),
        "markers": list(expression.columns),
        "directions": direction_metrics,
        "scope": "cross-severity and cross-context evidence; not human external validation",
    }
    (root / "cross_context_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True, sharey=True)
    panels = predictions.groupby(["model", "test_severity"], sort=True)
    for ax, ((model_name, test_severity), frame) in zip(axes.flat, panels):
        sns.scatterplot(data=frame, x="observed_hours", y="predicted_hours", s=90, ax=ax)
        ax.plot([0, 168], [0, 168], "--", color="#334155")
        ax.set_title(f"{model_name}: held-out {test_severity}")
        ax.set(xlabel="Observed time (hours)", ylabel="Predicted time (hours)")
    fig.suptitle("GSE162565 fixed-panel cross-severity validation")
    fig.tight_layout()
    fig.savefig(figures / "cross_severity_predictions.png", dpi=180)
    plt.close(fig)
    return summary
