from __future__ import annotations

import argparse
import math
import re
from collections import Counter

import numpy as np
import pandas as pd

from src.data.fetch_uniprot import UniProtRecord, fetch_uniprot_by_target_name, parse_fasta
from src.features.protein_grouping import charge_group, functional_group, hydrophobicity_group, pi_group
from src.utils.io import read_csv, write_csv
from src.utils.logging import get_logger

LOGGER = get_logger(__name__)

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:  # pragma: no cover
    ProteinAnalysis = None

ACIDIC = {"D", "E"}
BASIC = {"K", "R", "H"}
POLAR = {"S", "T", "N", "Q", "Y", "C", "W", "H", "K", "R", "D", "E"}
NONPOLAR = {"A", "V", "L", "I", "M", "F", "P", "G"}
AROMATIC = {"F", "W", "Y"}
AA_MASS = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
    "Q": 146.15, "E": 147.13, "G": 75.07, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}
KD_SCALE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def clean_protein_sequence(seq: object) -> str:
    if seq is None or pd.isna(seq):
        return ""
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", str(seq).upper())


def residue_fraction(counts: Counter[str], residues: set[str], length: int) -> float:
    return sum(counts.get(aa, 0) for aa in residues) / length if length else np.nan


def fallback_charge(seq: str) -> float:
    counts = Counter(seq)
    return float(counts.get("K", 0) + counts.get("R", 0) + 0.1 * counts.get("H", 0) - counts.get("D", 0) - counts.get("E", 0))


def compute_physical_features(seq: str) -> dict[str, float]:
    length = len(seq)
    counts = Counter(seq)
    if not seq:
        return {
            "protein_sequence_length": np.nan,
            "molecular_weight": np.nan,
            "theoretical_pI": np.nan,
            "charge_at_pH_7_4": np.nan,
            "GRAVY": np.nan,
            "acidic_residue_fraction": np.nan,
            "basic_residue_fraction": np.nan,
            "polar_residue_fraction": np.nan,
            "nonpolar_residue_fraction": np.nan,
            "aromatic_residue_fraction": np.nan,
            "cysteine_fraction": np.nan,
        }
    if ProteinAnalysis:
        analysis = ProteinAnalysis(seq)
        molecular_weight = analysis.molecular_weight()
        theoretical_pi = analysis.isoelectric_point()
        charge = analysis.charge_at_pH(7.4)
        gravy = analysis.gravy()
    else:
        molecular_weight = sum(AA_MASS.get(aa, 0.0) for aa in seq) - 18.015 * max(length - 1, 0)
        theoretical_pi = np.nan
        charge = fallback_charge(seq)
        gravy = sum(KD_SCALE.get(aa, 0.0) for aa in seq) / length
    return {
        "protein_sequence_length": length,
        "molecular_weight": molecular_weight,
        "theoretical_pI": theoretical_pi,
        "charge_at_pH_7_4": charge,
        "GRAVY": gravy,
        "acidic_residue_fraction": residue_fraction(counts, ACIDIC, length),
        "basic_residue_fraction": residue_fraction(counts, BASIC, length),
        "polar_residue_fraction": residue_fraction(counts, POLAR, length),
        "nonpolar_residue_fraction": residue_fraction(counts, NONPOLAR, length),
        "aromatic_residue_fraction": residue_fraction(counts, AROMATIC, length),
        "cysteine_fraction": counts.get("C", 0) / length,
    }


def lookup_sequence(target_name: str, fasta_records: dict[str, str]) -> tuple[str | None, str]:
    if target_name in fasta_records:
        return target_name, fasta_records[target_name]
    low = target_name.lower()
    for key, seq in fasta_records.items():
        if low == key.lower() or low in key.lower():
            return key, seq
    return None, ""


def build_protein_feature_table(
    pairs: pd.DataFrame,
    fasta_path: str | None = None,
    fetch_uniprot: bool = False,
    organism: str | None = "9606",
) -> pd.DataFrame:
    fasta_records = parse_fasta(fasta_path) if fasta_path else {}
    rows = []
    for target_name in sorted(pairs["target_name"].dropna().unique()):
        accession, seq = lookup_sequence(target_name, fasta_records)
        record = UniProtRecord(accession, target_name, seq)
        if not seq and fetch_uniprot:
            record = fetch_uniprot_by_target_name(target_name, organism=organism)
            seq = record.sequence or ""
        if not seq:
            LOGGER.warning("No protein sequence found for target: %s", target_name)
        annotation = " ".join([record.keywords or "", record.go_terms or "", record.comments or "", record.protein_name or ""])
        features = compute_physical_features(clean_protein_sequence(seq))
        rows.append(
            {
                "target_name": target_name,
                "uniprot_accession": record.accession,
                "uniprot_protein_name": record.protein_name,
                "protein_sequence": clean_protein_sequence(seq),
                **features,
                "charge_group": charge_group(features["charge_at_pH_7_4"]),
                "hydrophobicity_group": hydrophobicity_group(features["GRAVY"]),
                "pI_group": pi_group(features["theoretical_pI"]),
                "functional_group": functional_group(annotation),
                "annotation_text": annotation,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/aptamer_protein_pairs_clean.csv")
    parser.add_argument("--output", default="data/processed/protein_features.csv")
    parser.add_argument("--fasta", default=None, help="Optional local protein FASTA keyed by target name or accession.")
    parser.add_argument("--fetch-uniprot", action="store_true", help="Query UniProt when a sequence is not in FASTA.")
    parser.add_argument("--organism", default="9606", help="UniProt organism_id filter; use empty string to disable.")
    args = parser.parse_args()
    organism = args.organism or None
    pairs = read_csv(args.input)
    features = build_protein_feature_table(pairs, args.fasta, args.fetch_uniprot, organism)
    write_csv(features, args.output)
    LOGGER.info("Wrote %s protein rows to %s", len(features), args.output)


if __name__ == "__main__":
    main()
