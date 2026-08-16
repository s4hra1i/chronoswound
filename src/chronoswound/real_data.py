"""Reproducible exploratory analysis of the public GSE8056 dataset."""

from __future__ import annotations

import gzip
import json
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import f_oneway
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE8nnn/GSE8056/"
    "matrix/GSE8056_series_matrix.txt.gz"
)
ANNOTATION_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
)
GROUPS = ["0–3 days"] * 3 + ["4–7 days"] * 3 + [">7 days"] * 3 + ["Normal"] * 3
GROUP_ORDER = ["Normal", "0–3 days", "4–7 days", ">7 days"]
FOCUS_GENES = [
    "IL6", "TNF", "CXCL8", "MPO", "CD68", "VEGFA", "COL1A1", "MMP9", "TGFB1",
    "CCL2", "CXCL2", "PTGS2", "SERPINE1", "FN1", "KRT14",
]


def _download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        return
    request = urllib.request.Request(url, headers={"User-Agent": "ChronosWound/0.1"})
    with urllib.request.urlopen(request, timeout=90) as response, target.open("wb") as out:
        out.write(response.read())


def download_geo_files(data_dir: str | Path) -> tuple[Path, Path]:
    """Download processed expression and platform annotation from NCBI GEO."""
    root = Path(data_dir)
    matrix = root / "GSE8056_series_matrix.txt.gz"
    annotation = root / "GPL570.annot.gz"
    _download(MATRIX_URL, matrix)
    _download(ANNOTATION_URL, annotation)
    return matrix, annotation


def parse_series_matrix(path: str | Path) -> pd.DataFrame:
    """Read only the expression table from a GEO series-matrix file."""
    with gzip.open(path, "rt") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_begin"))
    end = next(i for i, line in enumerate(lines) if line.startswith("!series_matrix_table_end"))
    from io import StringIO

    table = pd.read_csv(StringIO("".join(lines[start + 1 : end])), sep="\t", quotechar='"')
    return table.set_index("ID_REF").apply(pd.to_numeric)


def parse_annotation(path: str | Path) -> pd.Series:
    """Return a probe-to-gene-symbol mapping, excluding ambiguous probes."""
    with gzip.open(path, "rt") as handle:
        lines = handle.readlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("!platform_table_begin"))
    from io import StringIO

    annotation = pd.read_csv(StringIO("".join(lines[start + 1 :])), sep="\t", low_memory=False)
    symbols = annotation.set_index("ID")["Gene symbol"].dropna().astype(str)
    symbols = symbols[~symbols.str.contains("///")]
    symbols = symbols[~symbols.isin({"---", "nan"})]
    return symbols


def collapse_to_genes(expression: pd.DataFrame, annotation: pd.Series) -> pd.DataFrame:
    """Log2-transform and collapse multiple probes by their median expression."""
    joined = expression.join(annotation.rename("gene"), how="inner")
    values = np.log2(joined.drop(columns="gene") + 1)
    values["gene"] = joined["gene"]
    return values.groupby("gene").median()


def _benjamini_hochberg(p_values: pd.Series) -> pd.Series:
    order = np.argsort(p_values.to_numpy())
    ranked = p_values.to_numpy()[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return pd.Series(result, index=p_values.index)


def differential_timecourse(genes: pd.DataFrame) -> pd.DataFrame:
    """Exploratory omnibus ANOVA across the four published sample groups."""
    group_index = pd.Series(GROUPS, index=genes.columns)
    p_values = genes.apply(
        lambda row: f_oneway(*(row[group_index == group].values for group in GROUP_ORDER)).pvalue,
        axis=1,
    )
    result = pd.DataFrame({"p_value": p_values, "fdr": _benjamini_hochberg(p_values)})
    means = {group: genes.loc[:, group_index == group].mean(axis=1) for group in GROUP_ORDER}
    for group, values in means.items():
        result[f"mean_{group}"] = values
    result["max_abs_log2fc_vs_normal"] = pd.concat(
        [(means[g] - means["Normal"]).abs() for g in GROUP_ORDER[1:]], axis=1
    ).max(axis=1)
    return result.sort_values(["fdr", "max_abs_log2fc_vs_normal"], ascending=[True, False])


def run_geo_analysis(
    matrix_path: str | Path, annotation_path: str | Path, output: str | Path
) -> dict:
    """Run QC, PCA, temporal testing and marker visualisation for GSE8056."""
    root = Path(output)
    figures = root / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    probes = parse_series_matrix(matrix_path)
    genes = collapse_to_genes(probes, parse_annotation(annotation_path))
    metadata = pd.DataFrame(
        {"sample": genes.columns, "time_bin": GROUPS, "replicate": [3, 1, 2] * 4}
    )
    metadata.to_csv(root / "sample_metadata.csv", index=False)
    stats = differential_timecourse(genes)
    stats.to_csv(root / "temporal_gene_statistics.csv")

    top_variable = genes.var(axis=1).nlargest(min(1000, len(genes))).index
    pca = PCA(n_components=2).fit_transform(
        StandardScaler().fit_transform(genes.loc[top_variable].T)
    )
    pca_df = metadata.assign(PC1=pca[:, 0], PC2=pca[:, 1])
    pca_df.to_csv(root / "pca_coordinates.csv", index=False)

    sns.set_theme(style="whitegrid", context="talk")
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.scatterplot(
        data=pca_df, x="PC1", y="PC2", hue="time_bin", hue_order=GROUP_ORDER,
        s=130, ax=ax, palette="viridis",
    )
    for _, row in pca_df.iterrows():
        ax.annotate(row["sample"], (row.PC1, row.PC2), fontsize=7, alpha=0.7)
    ax.set_title("GSE8056: PCA of 1,000 most variable genes")
    fig.tight_layout()
    fig.savefig(figures / "real_data_pca.png", dpi=180)
    plt.close(fig)

    present = [gene for gene in FOCUS_GENES if gene in genes.index]
    marker = genes.loc[present]
    marker_z = marker.sub(marker.mean(axis=1), axis=0).div(marker.std(axis=1), axis=0)
    labelled = marker_z.copy()
    labelled.columns = [f"{group}\nR{rep}" for group, rep in zip(GROUPS, [3, 1, 2] * 4)]
    grid = sns.clustermap(
        labelled, col_cluster=False, cmap="vlag", center=0, figsize=(12, 8),
        cbar_kws={"label": "Gene-wise z-score"},
    )
    grid.fig.suptitle("Literature-informed wound-response panel", y=1.02)
    grid.savefig(figures / "real_marker_heatmap.png", dpi=180)
    plt.close(grid.fig)

    top = stats.head(20).index
    top_z = genes.loc[top].sub(genes.loc[top].mean(axis=1), axis=0).div(
        genes.loc[top].std(axis=1), axis=0
    )
    grid = sns.clustermap(
        top_z, col_cluster=False, cmap="mako", figsize=(11, 9),
        xticklabels=[f"{g} R{r}" for g, r in zip(GROUPS, [3, 1, 2] * 4)],
    )
    grid.fig.suptitle("Top temporal signals (exploratory ANOVA)", y=1.02)
    grid.savefig(figures / "top_temporal_genes.png", dpi=180)
    plt.close(grid.fig)

    summary = {
        "accession": "GSE8056",
        "platform": "GPL570",
        "n_samples": int(genes.shape[1]),
        "n_genes": int(genes.shape[0]),
        "group_sizes": pd.Series(GROUPS).value_counts().to_dict(),
        "fdr_below_0.05": int((stats.fdr < 0.05).sum()),
        "analysis_status": "exploratory; pooled samples; no external validation",
    }
    (root / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    from .validation import evaluate_human_time_bins

    summary["fixed_panel_validation"] = evaluate_human_time_bins(genes, root)
    (root / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
