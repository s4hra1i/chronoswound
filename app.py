"""Interactive research explorer for all three ChronosWound evidence tracks."""

import json
from pathlib import Path

import joblib
import pandas as pd
import streamlit as st

from chronoswound.modelling import BASE_FEATURES

st.set_page_config(page_title="ChronosWound", page_icon="🧬", layout="wide")
st.title("ChronosWound research explorer")
st.warning("Research prototype only—not a clinical or forensic instrument.")

real_tab, cross_tab, demo_tab = st.tabs(
    ["Human transcriptomics", "Cross-context validation", "Synthetic estimator"]
)

with real_tab:
    human_summary = Path("reports/gse8056/analysis_summary.json")
    if human_summary.exists():
        human = json.loads(human_summary.read_text())
        validation = human["fixed_panel_validation"]
        a, b, c = st.columns(3)
        a.metric("Annotated genes", f"{human['n_genes']:,}")
        b.metric("Four-group LOOCV", f"{100 * validation['all_four_groups_accuracy']:.1f}%")
        c.metric(
            "Injured-only LOOCV",
            f"{100 * validation['injured_only_three_bins_accuracy']:.1f}%",
        )
        st.image("reports/gse8056/figures/real_data_pca.png")
        st.caption(
            "Twelve pooled arrays; confidence intervals are wide and results are not patient-level validation."
        )
    else:
        st.info("Run `chronoswound real-analysis` to populate this tab.")

with cross_tab:
    rat_summary = Path("reports/gse162565/cross_context_summary.json")
    if rat_summary.exists():
        rat = json.loads(rat_summary.read_text())
        directions = rat["directions"]
        left, right = st.columns(2)
        left.metric(
            "Ridge mild → severe MAE",
            f"{directions['train_mild_test_severe']['ridge']['mae_hours']:.1f} h",
        )
        right.metric(
            "Ridge severe → mild MAE",
            f"{directions['train_severe_test_mild']['ridge']['mae_hours']:.1f} h",
        )
        st.image("reports/gse162565/figures/cross_severity_predictions.png")
        st.caption("Rat muscle evidence tests severity robustness, not human transportability.")
    else:
        st.info("Run `chronoswound cross-context` to populate this tab.")

with demo_tab:
    model_path = Path("reports/model.joblib")
    if not model_path.exists():
        st.info("Run the synthetic training workflow first to create reports/model.joblib.")
    else:
        artifact = joblib.load(model_path)
        st.sidebar.header("Illustrative measurements")
        defaults = {name: 4.0 for name in BASE_FEATURES}
        defaults.update(
            {"neutrophil_pct": 35.0, "fibroblast_pct": 20.0, "temperature_c": 20.0}
        )
        values = {
            name: st.sidebar.number_input(name, value=defaults[name], step=0.1)
            for name in BASE_FEATURES
        }
        prediction = max(
            0.0, float(artifact["model"].predict(pd.DataFrame([values]))[0])
        )
        radius = float(artifact["interval_radius"])
        left, centre, right = st.columns(3)
        left.metric("Lower bound", f"{max(0, prediction - radius):.1f} h")
        centre.metric("Estimated age", f"{prediction:.1f} h")
        right.metric("Upper bound", f"{prediction + radius:.1f} h")
        st.caption(
            "The interval is calibrated on generated data and has no evidential meaning for real wounds."
        )
