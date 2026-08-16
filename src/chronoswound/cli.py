"""Command-line interface."""

from __future__ import annotations

import argparse
from pathlib import Path

from .data import generate_synthetic_cohort, load_dataset
from .modelling import train_and_evaluate
from .reporting import save_report
from .real_data import download_geo_files, run_geo_analysis
from .cross_context import download_rat_matrix, evaluate_cross_severity


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
    else:
        matrix = download_rat_matrix(args.data_dir)
        mapping = Path(__file__).parents[2] / "resources" / "GPL17117_focus_probe_map.csv"
        summary = evaluate_cross_severity(matrix, mapping, args.output)
        print(
            f"Cross-severity validation complete for "
            f"{summary['n_individual_wounded_animals']} wounded animals"
        )


if __name__ == "__main__":
    main()
