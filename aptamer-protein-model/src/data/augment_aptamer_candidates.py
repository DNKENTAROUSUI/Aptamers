from __future__ import annotations

import argparse

import pandas as pd

from src.data.load_aptamer_data import load_and_clean
from src.utils.io import write_csv
from src.utils.logging import get_logger

LOGGER = get_logger(__name__)

DNA_BASES = "ACGT"
RNA_BASES = "ACGU"


def next_base(base: str, alphabet: str, offset: int = 1) -> str:
    if base not in alphabet:
        return base
    return alphabet[(alphabet.index(base) + offset) % len(alphabet)]


def mutate_sequence(seq: str, variant_index: int, nucleic_acid_type: str) -> str:
    alphabet = RNA_BASES if nucleic_acid_type.upper() == "RNA" or "U" in seq else DNA_BASES
    if len(seq) < 8:
        return seq
    chars = list(seq)
    mutable_positions = list(range(3, len(chars) - 3))
    first = mutable_positions[(variant_index * 7) % len(mutable_positions)]
    second = mutable_positions[(variant_index * 11 + 5) % len(mutable_positions)]
    chars[first] = next_base(chars[first], alphabet, 1)
    if second != first:
        chars[second] = next_base(chars[second], alphabet, 2)
    return "".join(chars)


def build_candidate_variants(df: pd.DataFrame, variants_per_seed: int) -> pd.DataFrame:
    rows = []
    for _, row in df.iterrows():
        base = row.to_dict()
        base["candidate_status"] = base.get("candidate_status", "known_or_reported")
        base["parent_aptamer_id"] = base.get("parent_aptamer_id", "")
        rows.append(base)
        for idx in range(1, variants_per_seed + 1):
            variant = row.to_dict()
            variant["aptamer_id"] = f"{row['aptamer_id']}_VAR{idx:02d}"
            variant["aptamer_name"] = f"{row['aptamer_name']}-candidate-v{idx}"
            variant["sequence"] = mutate_sequence(str(row["sequence"]), idx, str(row["nucleic_acid_type"]))
            variant["kd_value"] = ""
            variant["kd_unit"] = ""
            variant["source_database"] = "exploratory_variant"
            variant["reference_doi"] = ""
            variant["candidate_status"] = "in_silico_variant"
            variant["parent_aptamer_id"] = row["aptamer_id"]
            rows.append(variant)
    out = pd.DataFrame(rows)
    out = out.drop_duplicates(subset=["sequence", "target_name", "candidate_status"]).reset_index(drop=True)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Create deterministic exploratory aptamer candidate variants.")
    parser.add_argument("--input", default="data/raw/aptamers_literature_seed.csv")
    parser.add_argument("--output", default="data/raw/aptamers_candidate_expanded.csv")
    parser.add_argument("--variants-per-seed", type=int, default=2)
    args = parser.parse_args()
    if args.variants_per_seed < 0:
        raise ValueError("--variants-per-seed must be non-negative")
    seed = load_and_clean(args.input)
    expanded = build_candidate_variants(seed, args.variants_per_seed)
    write_csv(expanded, args.output)
    LOGGER.info("Wrote %s rows to %s", len(expanded), args.output)


if __name__ == "__main__":
    main()
