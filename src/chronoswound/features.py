"""Biologically motivated, transparent feature engineering."""

from __future__ import annotations

from typing import ClassVar

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin


class PhaseScoreTransformer(BaseEstimator, TransformerMixin):
    """Add mean expression scores for broad wound-healing phases."""

    phase_genes: ClassVar[dict[str, list[str]]] = {
        "inflammation_score": ["IL6", "TNF", "CXCL8", "MPO"],
        "repair_score": ["CD68", "VEGFA", "MMP9"],
        "remodelling_score": ["COL1A1", "TGFB1"],
    }

    def fit(self, X: pd.DataFrame, y=None):
        self.feature_names_in_ = np.asarray(X.columns, dtype=object)
        return self

    def transform(self, X: pd.DataFrame) -> pd.DataFrame:
        out = X.copy()
        for score, genes in self.phase_genes.items():
            out[score] = out[genes].mean(axis=1)
        out["immune_balance"] = out["neutrophil_pct"] - out["fibroblast_pct"]
        return out

    def get_feature_names_out(self, input_features=None):
        """Expose generated columns to downstream scikit-learn components."""
        if input_features is None:
            input_features = self.feature_names_in_
        added = [*self.phase_genes, "immune_balance"]
        return np.asarray([*input_features, *added], dtype=object)
