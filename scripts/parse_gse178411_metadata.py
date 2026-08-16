"""Parse GSE178411 sample characteristics from a GEO SOFT family file."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def parse_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" else open
    samples: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    wound_type_seen = 0

    with opener(path, "rt", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {"geo_accession": line.split(" = ", 1)[1]}
                wound_type_seen = 0
            elif current is None:
                continue
            elif line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_characteristics_ch1 = "):
                characteristic = line.split(" = ", 1)[1]
                key, value = characteristic.split(": ", 1)
                key = key.strip().lower().replace(" ", "_")
                if key == "wound_type":
                    key = "wound_stage" if wound_type_seen == 0 else "sample_class"
                    wound_type_seen += 1
                current[key] = value.strip()

    if current is not None:
        samples.append(current)

    metadata = pd.DataFrame(samples)
    metadata["patient_id"] = metadata["subject"].str.extract(r"(\d+)").astype("Int64")
    metadata["days_since_injury"] = pd.to_numeric(
        metadata["days_since_injury"].replace("--", pd.NA), errors="coerce"
    ).astype("Int64")
    metadata["age"] = pd.to_numeric(metadata["age"], errors="coerce").astype("Float64")
    return metadata[
        [
            "geo_accession",
            "title",
            "patient_id",
            "wound_stage",
            "sample_class",
            "days_since_injury",
            "age",
            "sex",
            "hispanic",
            "race",
            "burn_type",
            "location",
            "tissue",
        ]
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("soft", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    metadata = parse_soft(args.soft)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.output, index=False)
    print(f"Wrote {len(metadata)} samples to {args.output}")


if __name__ == "__main__":
    main()
