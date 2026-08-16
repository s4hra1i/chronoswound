"""Command-line interface."""

from __future__ import annotations

import argparse
from pathlib import Path

from .cross_context import download_rat_matrix, evaluate_cross_severity
from .data import generate_synthetic_cohort, load_dataset
from .gse178411 import AnalysisConfig, download_counts, run_analysis
from .modelling import train_and_evaluate
from .real_data import download_geo_files, run_geo_analysis
from .reporting import save_report


def main() -> None:
    parser = argparse.ArgumentParser(prog="chronoswound")
    sub = parser.add_subparsers(dest="command", required=True)
    generate = sub.add_parser("generate", help="Generate a synthetic demonstration cohort")
    generate.add_argument("--samples", type=int, default=360)
    generate.add_argument("--seed", type=int, default=42)
    generate.add_argument("--output", type=Path, default=Path("data/synthetic_wounds.csv"))
    train = sub.add_parser("train", help="Train, evaluate and report a wound-age model")
    train.add_argument("--input", type=Path, required=True)
    train.add_argument("--output", type=Path, default=Path("reports"))
    train.add_argument("--seed", type=int, default=42)
    real = sub.add_parser("real-analysis", help="Analyse the public GSE8056 burn dataset")
    real.add_argument("--data-dir", type=Path, default=Path("data/geo"))
    real.add_argument("--output", type=Path, default=Path("reports/gse8056"))
    cross = sub.add_parser(
        "cross-context", help="Run GSE162565 rat cross-severity validation"
    )
    cross.add_argument("--data-dir", type=Path, default=Path("data/geo"))
    cross.add_argument("--output", type=Path, default=Path("reports/gse162565"))
    prospective = sub.add_parser(
        "gse178411", help="Run the prospectively tagged GSE178411 analysis"
    )
    prospective.add_argument(
        "--counts", type=Path, default=Path("data/geo/GSE178411_counts.txt.gz")
    )
    prospective.add_argument(
        "--metadata", type=Path, default=Path("data/GSE178411_sample_metadata.csv")
    )
    prospective.add_argument("--output", type=Path, default=Path("reports/gse178411"))
    prospective.add_argument("--bootstrap-replicates", type=int, default=10_000)
    prospective.add_argument("--permutation-replicates", type=int, default=199)
    args = parser.parse_args()

    if args.command == "generate":
        args.output.parent.mkdir(parents=True, exist_ok=True)
        generate_synthetic_cohort(args.samples, args.seed).to_csv(args.output, index=False)
        print(f"Wrote {args.samples} synthetic observations to {args.output}")
    elif args.command == "train":
        result = train_and_evaluate(load_dataset(args.input), args.seed)
        save_report(result, args.output)
        print(f"Selected {result.metrics['selected_model']}; test MAE: {result.metrics['test_mae_hours']:.1f} hours")
        print(f"Report written to {args.output}")
    elif args.command == "real-analysis":
        matrix, annotation = download_geo_files(args.data_dir)
        summary = run_geo_analysis(matrix, annotation, args.output)
        print(
            f"Analysed GSE8056: {summary['n_samples']} arrays and "
            f"{summary['n_genes']} annotated genes"
        )
        print(f"Exploratory report written to {args.output}")
    elif args.command == "cross-context":
        matrix = download_rat_matrix(args.data_dir)
        mapping = Path(__file__).parents[2] / "resources" / "GPL17117_focus_probe_map.csv"
        summary = evaluate_cross_severity(matrix, mapping, args.output)
        print(
            f"Cross-severity validation complete for "
            f"{summary['n_individual_wounded_animals']} wounded animals"
        )
    else:
        counts = download_counts(args.counts)
        config = AnalysisConfig(
            bootstrap_replicates=args.bootstrap_replicates,
            permutation_replicates=args.permutation_replicates,
        )
        summary = run_analysis(counts, args.metadata, args.output, config)
        result = summary["primary_comparison"]
        print(
            "GSE178411 analysis complete; prospective success rule met: "
            f"{result['all_three_success_conditions_met']}"
        )
        print(f"Report written to {args.output}")


if __name__ == "__main__":
    main()
