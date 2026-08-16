"""Prospectively specified GSE178411 early/late-subacute wound analysis."""

from __future__ import annotations

import json
import os
import platform
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import sklearn
from joblib import Parallel, delayed
from scipy.stats import norm
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import accuracy_score, balanced_accuracy_score, mean_absolute_error
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

PROTOCOL_TAG = "gse178411-protocol-v1.0"
PROTOCOL_SEED = 20260816
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178411/suppl/"
    "GSE178411_counts.txt.gz"
)
PANEL_GENE_IDS = {
    "IL6": 3569,
    "TNF": 7124,
    "CXCL8": 3576,
    "CXCL2": 2920,
    "PTGS2": 5743,
    "MPO": 4353,
    "CD68": 968,
    "CCL2": 6347,
    "MMP9": 4318,
    "VEGFA": 7422,
    "TGFB1": 7040,
    "COL1A1": 1277,
    "FN1": 2335,
    "SERPINE1": 5054,
    "KRT14": 3861,
}
PANEL = list(PANEL_GENE_IDS)
COVARIATES = ["age", "sex", "burn_type", "location"]
CATEGORICAL = ["sex", "burn_type", "location"]
ALPHAS = np.asarray([10 ** (-3 + 0.25 * k) for k in range(25)])


@dataclass(frozen=True)
class AnalysisConfig:
    """Runtime settings whose defaults reproduce protocol version 1.0."""

    seeds: tuple[int, ...] = tuple(range(PROTOCOL_SEED, PROTOCOL_SEED + 20))
    outer_splits: int = 5
    inner_splits: int = 5
    bootstrap_replicates: int = 10_000
    permutation_replicates: int = 199
    permutation_jobs: int = -1


def download_counts(destination: Path) -> Path:
    """Download the study-author-submitted filtered count matrix if absent."""
    if not destination.exists():
        destination.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(COUNTS_URL, destination)
    return destination


def load_counts(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    """Load submitted counts and return locked genes plus sample-local library sizes."""
    opener = "gzip" if path.suffix == ".gz" else "infer"
    with pd.io.common.get_handle(path, "r", compression=opener) as handle:
        labels = handle.handle.readline().rstrip("\n").split("\t")
    raw = pd.read_csv(path, sep="\t", header=None, skiprows=1, compression=opener)
    if raw.shape[1] != len(labels) + 1:
        raise ValueError("Unexpected GSE178411 count-matrix shape")
    gene_ids = pd.to_numeric(raw.iloc[:, 0], errors="raise").astype(int)
    values = raw.iloc[:, 1:].astype(float)
    values.columns = labels
    library_sizes = values.sum(axis=0)
    if (library_sizes <= 0).any():
        raise ValueError("Count matrix contains an empty library")
    values.index = gene_ids
    missing = sorted(set(PANEL_GENE_IDS.values()) - set(values.index))
    if missing:
        raise ValueError(f"Locked Entrez IDs missing from counts: {missing}")
    panel = values.loc[list(PANEL_GENE_IDS.values())].copy()
    panel.index = PANEL
    return panel.T, library_sizes


def prepare_cohort(counts_path: Path, metadata_path: Path) -> pd.DataFrame:
    """Apply the prospectively fixed eligibility rules and sample-local log-CPM."""
    panel_counts, library_sizes = load_counts(counts_path)
    expression = np.log2(panel_counts.div(library_sizes, axis=0) * 1_000_000 + 0.5)
    expression.index.name = "sample_label"
    expression = expression.reset_index()

    metadata = pd.read_csv(metadata_path)
    metadata["sample_label"] = metadata["title"].str.split(":", n=1).str[0]
    eligible = metadata.loc[
        metadata["sample_class"].eq("wound")
        & metadata["wound_stage"].isin(["Early Wound", "Late wound"])
        & metadata["days_since_injury"].notna()
        & metadata["age"].notna()
    ].copy()
    cohort = eligible.merge(expression, on="sample_label", how="left", validate="one_to_one")
    if cohort[PANEL].isna().any().any():
        raise ValueError("Eligible metadata samples did not all match the count matrix")
    if len(cohort) != 49 or cohort["patient_id"].nunique() != 39:
        raise ValueError("Prospective complete-case cohort must contain 49 samples/39 patients")
    cohort["patient_id"] = cohort["patient_id"].astype(int)
    cohort["days_since_injury"] = cohort["days_since_injury"].astype(float)
    return cohort.sort_values(["patient_id", "geo_accession"]).reset_index(drop=True)


def _preprocessor(model: str) -> ColumnTransformer:
    transformers = []
    if model in {"covariates", "combined"}:
        numeric = Pipeline(
            [("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]
        )
        categorical = Pipeline(
            [
                ("impute", SimpleImputer(strategy="most_frequent")),
                ("encode", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
                ("scale", StandardScaler()),
            ]
        )
        transformers.extend(
            [("age", numeric, ["age"]), ("categorical", categorical, CATEGORICAL)]
        )
    if model in {"panel", "combined"}:
        transformers.append(("panel", StandardScaler(), PANEL))
    return ColumnTransformer(transformers)


def _fit_predict(
    model: str,
    train: pd.DataFrame,
    test: pd.DataFrame,
    inner_seed: int,
    inner_splits: int,
) -> tuple[np.ndarray, float]:
    y_train = train["days_since_injury"].to_numpy()
    groups = train["patient_id"].to_numpy()
    cv = GroupKFold(n_splits=inner_splits, shuffle=True, random_state=inner_seed)
    errors = np.zeros(len(ALPHAS))
    for inner_train, inner_test in cv.split(train, groups=groups):
        preprocessing = _preprocessor(model)
        x_train = preprocessing.fit_transform(train.iloc[inner_train])
        x_test = preprocessing.transform(train.iloc[inner_test])
        for alpha_index, alpha in enumerate(ALPHAS):
            fitted = Ridge(alpha=alpha).fit(x_train, y_train[inner_train])
            errors[alpha_index] += mean_absolute_error(
                y_train[inner_test], fitted.predict(x_test)
            )
    best_alpha = float(ALPHAS[np.argmin(errors)])
    final = Pipeline(
        [("preprocess", _preprocessor(model)), ("ridge", Ridge(alpha=best_alpha))]
    )
    final.fit(train, y_train)
    return final.predict(test), best_alpha


def repeated_nested_cv(
    cohort: pd.DataFrame,
    config: AnalysisConfig,
    outcome: np.ndarray | None = None,
    retain_rows: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run the four-model repeated nested grouped evaluation."""
    work = cohort.copy()
    if outcome is not None:
        work["days_since_injury"] = outcome
    predictions: list[dict[str, object]] = []
    folds: list[dict[str, object]] = []
    assignments: list[dict[str, object]] = []
    groups = work["patient_id"].to_numpy()

    for repetition, seed in enumerate(config.seeds):
        outer = GroupKFold(
            n_splits=config.outer_splits, shuffle=True, random_state=seed
        )
        for fold, (train_idx, test_idx) in enumerate(outer.split(work, groups=groups)):
            train = work.iloc[train_idx]
            test = work.iloc[test_idx]
            observed = test["days_since_injury"].to_numpy()
            fold_predictions: dict[str, np.ndarray] = {
                "median": np.repeat(train["days_since_injury"].median(), len(test))
            }
            selected_alphas: dict[str, float | None] = {"median": None}
            for model in ("covariates", "panel", "combined"):
                pred, alpha = _fit_predict(
                    model,
                    train,
                    test,
                    inner_seed=seed + 10_000 + fold,
                    inner_splits=config.inner_splits,
                )
                fold_predictions[model] = pred
                selected_alphas[model] = alpha

            for model, predicted in fold_predictions.items():
                folds.append(
                    {
                        "repetition": repetition,
                        "seed": seed,
                        "fold": fold,
                        "model": model,
                        "n_test": len(test),
                        "mae_days": mean_absolute_error(observed, predicted),
                        "selected_alpha": selected_alphas[model],
                    }
                )
                if retain_rows:
                    for row_number, (_, row) in enumerate(test.iterrows()):
                        predictions.append(
                            {
                                "repetition": repetition,
                                "seed": seed,
                                "fold": fold,
                                "model": model,
                                "geo_accession": row["geo_accession"],
                                "patient_id": int(row["patient_id"]),
                                "observed_days": float(observed[row_number]),
                                "predicted_days": float(predicted[row_number]),
                                "absolute_error_days": float(
                                    abs(observed[row_number] - predicted[row_number])
                                ),
                            }
                        )
            if retain_rows:
                for patient in sorted(test["patient_id"].unique()):
                    assignments.append(
                        {
                            "repetition": repetition,
                            "seed": seed,
                            "patient_id": int(patient),
                            "fold": fold,
                        }
                    )
    return pd.DataFrame(predictions), pd.DataFrame(folds), pd.DataFrame(assignments)


def _repetition_metrics(predictions: pd.DataFrame) -> pd.DataFrame:
    return (
        predictions.groupby(["repetition", "seed", "model"], as_index=False)[
            "absolute_error_days"
        ]
        .mean()
        .rename(columns={"absolute_error_days": "pooled_mae_days"})
    )


def _paired_patient_values(
    predictions: pd.DataFrame, reference: str, candidate: str
) -> pd.DataFrame:
    keys = ["repetition", "geo_accession", "patient_id"]
    wide = predictions.pivot(index=keys, columns="model", values="absolute_error_days")
    wide["difference"] = wide[reference] - wide[candidate]
    return wide.reset_index()


def bca_cluster_interval(
    paired: pd.DataFrame,
    replicates: int,
    seed: int,
) -> dict[str, float | int]:
    """BCa interval for mean paired error difference, resampling patients."""
    patient_values = {
        patient: group["difference"].to_numpy()
        for patient, group in paired.groupby("patient_id")
    }
    patients = np.asarray(sorted(patient_values))
    observed = float(paired["difference"].mean())
    rng = np.random.default_rng(seed)
    boot = np.empty(replicates)
    for index in range(replicates):
        selected = rng.choice(patients, size=len(patients), replace=True)
        boot[index] = np.concatenate([patient_values[p] for p in selected]).mean()
    proportion = np.clip(np.mean(boot < observed), 1 / (2 * replicates), 1 - 1 / (2 * replicates))
    z0 = norm.ppf(proportion)
    jack = np.asarray(
        [paired.loc[paired["patient_id"] != patient, "difference"].mean() for patient in patients]
    )
    jack_mean = jack.mean()
    numerator = np.sum((jack_mean - jack) ** 3)
    denominator = 6 * np.sum((jack_mean - jack) ** 2) ** 1.5
    acceleration = float(numerator / denominator) if denominator else 0.0

    adjusted = []
    for probability in (0.025, 0.975):
        z = norm.ppf(probability)
        adjusted.append(norm.cdf(z0 + (z0 + z) / (1 - acceleration * (z0 + z))))
    lower, upper = np.quantile(boot, np.clip(adjusted, 0, 1))
    return {
        "estimate_days": observed,
        "lower_95_days": float(lower),
        "upper_95_days": float(upper),
        "replicates": replicates,
        "patients": len(patients),
    }


def _permuted_outcome(cohort: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    outcome = cohort["days_since_injury"].to_numpy().copy()
    patient_rows = {
        patient: np.asarray(indices)
        for patient, indices in cohort.groupby("patient_id", sort=True).indices.items()
    }
    strata: dict[int, list[int]] = {}
    for patient, indices in patient_rows.items():
        strata.setdefault(len(indices), []).append(patient)
    for patients in strata.values():
        donors = rng.permutation(patients)
        original = outcome.copy()
        for recipient, donor in zip(patients, donors, strict=True):
            outcome[patient_rows[recipient]] = original[patient_rows[donor]]
    return outcome


def permutation_test(
    cohort: pd.DataFrame,
    config: AnalysisConfig,
    observed_improvement: float,
) -> tuple[dict[str, float | int], pd.DataFrame]:
    rng = np.random.default_rng(PROTOCOL_SEED + 900_000)
    outcomes = [
        _permuted_outcome(cohort, rng) for _ in range(config.permutation_replicates)
    ]

    def evaluate(permutation: int, outcome: np.ndarray) -> dict[str, float | int]:
        _, folds, _ = repeated_nested_cv(
            cohort, config, outcome=outcome, retain_rows=False
        )
        folds = folds.assign(weighted_error=folds["mae_days"] * folds["n_test"])
        repetition = folds.groupby(["repetition", "model"], as_index=False).agg(
            weighted_error=("weighted_error", "sum"), n_test=("n_test", "sum")
        )
        repetition["pooled_mae_days"] = repetition["weighted_error"] / repetition["n_test"]
        means = repetition.groupby("model")["pooled_mae_days"].mean()
        return {
            "permutation": permutation,
            "baseline_minus_panel_mae_days": float(means["median"] - means["panel"]),
        }

    rows = Parallel(n_jobs=config.permutation_jobs, verbose=10)(
        delayed(evaluate)(index, outcome) for index, outcome in enumerate(outcomes)
    )
    null = pd.DataFrame(
        rows, columns=["permutation", "baseline_minus_panel_mae_days"]
    )
    exceedances = int(
        (null["baseline_minus_panel_mae_days"] >= observed_improvement).sum()
    )
    return (
        {
            "p_value_one_sided": (1 + exceedances) / (1 + config.permutation_replicates),
            "replicates": config.permutation_replicates,
            "exceedances": exceedances,
        },
        null,
    )


def _distribution(values: pd.Series) -> dict[str, float | list[float]]:
    return {
        "mean": float(values.mean()),
        "median": float(values.median()),
        "iqr": [float(values.quantile(0.25)), float(values.quantile(0.75))],
        "range": [float(values.min()), float(values.max())],
        "values": [float(value) for value in values],
    }


def _save_figure(repetition_metrics: pd.DataFrame, output: Path) -> None:
    order = ["median", "covariates", "panel", "combined"]
    data = [
        repetition_metrics.loc[
            repetition_metrics["model"].eq(model), "pooled_mae_days"
        ].to_numpy()
        for model in order
    ]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot(data, tick_labels=["Median", "Covariates", "Panel", "Combined"])
    ax.set_ylabel("Pooled out-of-fold MAE (days)")
    ax.set_title("GSE178411 repeated patient-grouped validation")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def run_analysis(
    counts_path: Path,
    metadata_path: Path,
    output_dir: Path,
    config: AnalysisConfig | None = None,
) -> dict[str, object]:
    """Run and persist the complete prospective analysis."""
    if config is None:
        config = AnalysisConfig()
    cohort = prepare_cohort(counts_path, metadata_path)
    predictions, fold_metrics, assignments = repeated_nested_cv(cohort, config)
    repetition_metrics = _repetition_metrics(predictions)
    model_means = repetition_metrics.groupby("model")["pooled_mae_days"].mean()

    baseline_panel = _paired_patient_values(predictions, "median", "panel")
    covariates_combined = _paired_patient_values(predictions, "covariates", "combined")
    panel_interval = bca_cluster_interval(
        baseline_panel, config.bootstrap_replicates, PROTOCOL_SEED + 700_000
    )
    incremental_interval = bca_cluster_interval(
        covariates_combined, config.bootstrap_replicates, PROTOCOL_SEED + 800_000
    )
    observed_improvement = float(model_means["median"] - model_means["panel"])
    permutation, null = permutation_test(cohort, config, observed_improvement)
    relative_reduction = observed_improvement / float(model_means["median"])

    panel_success = bool(
        relative_reduction >= 0.20
        and panel_interval["lower_95_days"] > 0
        and permutation["p_value_one_sided"] < 0.05
    )
    incremental_success = bool(incremental_interval["lower_95_days"] > 0)

    panel_rows = predictions.loc[predictions["model"].eq("panel")].copy()
    panel_rows["observed_class"] = np.where(panel_rows["observed_days"] < 8, "Early", "Late")
    panel_rows["predicted_class"] = np.where(panel_rows["predicted_days"] < 8, "Early", "Late")

    summary: dict[str, object] = {
        "protocol_tag": PROTOCOL_TAG,
        "generated_from_commit": os.environ.get("CHRONOSWOUND_COMMIT_SHA", "not-recorded"),
        "cohort": {
            "samples": len(cohort),
            "patients": cohort["patient_id"].nunique(),
            "day_range": [
                float(cohort["days_since_injury"].min()),
                float(cohort["days_since_injury"].max()),
            ],
            "excluded_missing_age": ["GSM5390619"],
        },
        "configuration": {
            **asdict(config),
            "seeds": list(config.seeds),
            "alpha_grid": [float(alpha) for alpha in ALPHAS],
            "target": "untransformed days",
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "model_mae_days": {
            model: _distribution(
                repetition_metrics.loc[
                    repetition_metrics["model"].eq(model), "pooled_mae_days"
                ]
            )
            for model in ("median", "covariates", "panel", "combined")
        },
        "primary_comparison": {
            "baseline_minus_panel_mae": panel_interval,
            "relative_mae_reduction": float(relative_reduction),
            "permutation": permutation,
            "all_three_success_conditions_met": panel_success,
        },
        "incremental_comparison": {
            "covariates_minus_combined_mae": incremental_interval,
            "incremental_molecular_value_established": incremental_success,
        },
        "secondary_panel_threshold_classification": {
            "qualification": "Descriptive only; labels are a deterministic 7/8-day threshold.",
            "accuracy": accuracy_score(
                panel_rows["observed_class"], panel_rows["predicted_class"]
            ),
            "balanced_accuracy": balanced_accuracy_score(
                panel_rows["observed_class"], panel_rows["predicted_class"]
            ),
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(output_dir / "oof_predictions.csv", index=False)
    fold_metrics.to_csv(output_dir / "fold_metrics.csv", index=False)
    assignments.to_csv(output_dir / "fold_assignments.csv", index=False)
    repetition_metrics.to_csv(output_dir / "repetition_metrics.csv", index=False)
    null.to_csv(output_dir / "permutation_null.csv", index=False)
    cohort[["geo_accession", "patient_id", "days_since_injury"]].to_csv(
        output_dir / "analysis_cohort.csv", index=False
    )
    with (output_dir / "analysis_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    _save_figure(repetition_metrics, output_dir / "figures" / "model_mae_comparison.png")
    return summary
