from __future__ import annotations

import argparse
import re

import pandas as pd

from src.data.normalize_targets import normalize_target_name, normalize_target_type
from src.utils.io import missing_columns, read_csv, write_csv
from src.utils.logging import get_logger

LOGGER = get_logger(__name__)

REQUIRED_COLUMNS = [
    "aptamer_id",
    "aptamer_name",
    "sequence",
    "nucleic_acid_type",
    "target_name",
    "target_type",
    "kd_value",
    "kd_unit",
    "source_database",
    "reference_doi",
]


def normalize_sequence(value: object) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"[^ACGTUacgtu]", "", str(value)).upper()


def load_and_clean(input_path: str) -> pd.DataFrame:
    df = read_csv(input_path)
    missing = missing_columns(df, REQUIRED_COLUMNS)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    clean = df.copy()
    clean["sequence"] = clean["sequence"].map(normalize_sequence)
    clean["target_name"] = clean["target_name"].map(normalize_target_name)
    clean["target_type"] = clean["target_type"].map(normalize_target_type)
    clean["nucleic_acid_type"] = clean["nucleic_acid_type"].fillna("").astype(str).str.upper()
    clean["kd_value"] = pd.to_numeric(clean["kd_value"], errors="coerce")
    before = len(clean)
    clean = clean[(clean["target_type"] == "protein") & clean["sequence"].ne("") & clean["target_name"].ne("")]
    clean = clean.drop_duplicates(subset=["aptamer_id", "sequence", "target_name"]).reset_index(drop=True)
    LOGGER.info("Loaded %s rows; kept %s protein aptamer-target pairs", before, len(clean))
    return clean


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/raw/aptamers.csv")
    parser.add_argument("--output", default="data/processed/aptamer_protein_pairs_clean.csv")
    args = parser.parse_args()
    write_csv(load_and_clean(args.input), args.output)
    LOGGER.info("Wrote cleaned data to %s", args.output)


if __name__ == "__main__":
    main()
