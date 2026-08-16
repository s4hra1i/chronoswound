"""Dataset generation and validation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

GENES = ["IL6", "TNF", "CXCL8", "MPO", "CD68", "VEGFA", "COL1A1", "MMP9", "TGFB1"]
REQUIRED_COLUMNS = ["donor_id", "wound_age_hours", *GENES, "neutrophil_pct", "fibroblast_pct", "temperature_c"]


def _peak(age: np.ndarray, centre: float, width: float, height: float) -> np.ndarray:
    return height * np.exp(-0.5 * ((age - centre) / width) ** 2)


def generate_synthetic_cohort(n_samples: int = 360, seed: int = 42) -> pd.DataFrame:
    """Create a demonstrator cohort with donor and assay variation.

    Temporal curves are illustrative, not empirical reference ranges.
    """
    if n_samples < 60:
        raise ValueError("At least 60 samples are required for grouped modelling.")
    rng = np.random.default_rng(seed)
    n_donors = max(30, n_samples // 4)
    donor = rng.integers(1, n_donors + 1, size=n_samples)
    age = np.clip(rng.gamma(shape=1.7, scale=42, size=n_samples), 0.5, 240)
    donor_effects = rng.normal(0, 0.35, n_donors + 1)[donor]
    noise = lambda scale=0.45: rng.normal(0, scale, n_samples)

    data = {
        "donor_id": [f"D{x:03d}" for x in donor],
        "wound_age_hours": age.round(2),
        "IL6": 3.0 + _peak(age, 8, 10, 5.2) + donor_effects + noise(),
        "TNF": 2.8 + _peak(age, 12, 14, 4.4) + donor_effects + noise(),
        "CXCL8": 3.2 + _peak(age, 18, 18, 5.0) + donor_effects + noise(),
        "MPO": 2.5 + _peak(age, 24, 22, 4.7) + donor_effects + noise(),
        "CD68": 2.4 + _peak(age, 54, 35, 4.0) + donor_effects + noise(),
        "VEGFA": 2.5 + _peak(age, 76, 45, 3.7) + donor_effects + noise(),
        "COL1A1": 2.0 + 4.5 / (1 + np.exp(-(age - 92) / 24)) + donor_effects + noise(),
        "MMP9": 2.5 + _peak(age, 48, 32, 4.2) + donor_effects + noise(),
        "TGFB1": 2.2 + 3.8 / (1 + np.exp(-(age - 64) / 22)) + donor_effects + noise(),
        "neutrophil_pct": np.clip(8 + _peak(age, 20, 24, 66) + noise(5), 0, 100),
        "fibroblast_pct": np.clip(5 + 50 / (1 + np.exp(-(age - 80) / 22)) + noise(4), 0, 100),
        "temperature_c": rng.normal(20, 3.5, n_samples),
    }
    return pd.DataFrame(data).round(3)


def validate_dataset(df: pd.DataFrame) -> None:
    missing = sorted(set(REQUIRED_COLUMNS) - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    if df[REQUIRED_COLUMNS].isna().any().any():
        raise ValueError("Required columns contain missing values.")
    if (df["wound_age_hours"] < 0).any():
        raise ValueError("Wound age cannot be negative.")
    if df["donor_id"].nunique() < 10:
        raise ValueError("At least 10 unique donors are required.")


def load_dataset(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    validate_dataset(df)
    return df
