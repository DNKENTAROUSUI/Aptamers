from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from src.data.load_aptamer_data import REQUIRED_COLUMNS, load_and_clean
from src.utils.io import write_csv
from src.utils.logging import get_logger

LOGGER = get_logger(__name__)


def merge_sources(inputs: list[str]) -> pd.DataFrame:
    frames = []
    for path in inputs:
        source = Path(path)
        if not source.exists():
            LOGGER.warning("Skipping missing source: %s", path)
            continue
        frames.append(load_and_clean(str(source)))
    if not frames:
        return pd.DataFrame(columns=REQUIRED_COLUMNS)
    merged = pd.concat(frames, ignore_index=True)
    merged = merged.drop_duplicates(subset=["sequence", "target_name"]).reset_index(drop=True)
    return merged


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge multiple aptamer CSV sources into one cleaned table.")
    parser.add_argument("--inputs", nargs="+", required=True)
    parser.add_argument("--output", default="data/processed/aptamer_protein_pairs_clean.csv")
    args = parser.parse_args()
    merged = merge_sources(args.inputs)
    write_csv(merged, args.output)
    LOGGER.info("Merged %s rows from %s sources into %s", len(merged), len(args.inputs), args.output)


if __name__ == "__main__":
    main()
