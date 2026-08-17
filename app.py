"""Interactive research explorer for the four ChronosWound evidence tracks."""

import json
from pathlib import Path

import joblib
import pandas as pd
import streamlit as st

from chronoswound.modelling import BASE_FEATURES

st.set_page_config(page_title="ChronosWound", page_icon="🧬", layout="wide")
st.title("ChronosWound research explorer")
st.warning("Research prototype only—not a clinical or forensic instrument.")

primary_tab, real_tab, cross_tab, demo_tab = st.tabs(
    [
        "Primary human validation",
        "Exploratory human arrays",
        "Cross-context validation",
        "Synthetic demonstration",
    ]
)

with primary_tab:
    primary_summary = Path("reports/gse178411/analysis_summary.json")
    primary_figure = Path("reports/gse178411/figures/model_mae_comparison.png")
    if primary_summary.exists():
        primary = json.loads(primary_summary.read_text())
        model_metrics = primary["model_mae_days"]
        comparison = primary["primary_comparison"]
        a, b, c = st.columns(3)
        a.metric("Locked-panel MAE", f"{model_metrics['panel']['mean']:.2f} days")
        b.metric("Training-median MAE", f"{model_metrics['median']['mean']:.2f} days")
        c.metric(
            "Relative MAE reduction",
            f"{100 * comparison['relative_mae_reduction']:.1f}%",
        )
        interval = comparison["baseline_minus_panel_mae"]
        st.write(
            "Patient-clustered improvement: "
            f"{interval['estimate_days']:.2f} days "
            f"(95% CI {interval['lower_95_days']:.2f}–"
            f"{interval['upper_95_days']:.2f})."
        )
        if primary_figure.exists():
            st.image(str(primary_figure))
        st.caption(
            "Internal validation in 49 samples from 39 patients. "
            "This is not external, clinical or forensic validation."
        )
    else:
        st.info("Run `chronoswound gse178411` to populate this tab.")

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
    model_path = Path("reports/synthetic/model.joblib")
    if not model_path.exists():
        st.info(
            "Run the synthetic training workflow first to create "
            "reports/synthetic/model.joblib."
        )
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
